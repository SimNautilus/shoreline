// std includes
#include <iostream>
#include <string>
#include <vector>
#include <chrono>

// mpi includes
#include <mpi.h>

// ShorelineDDG includes
#include "MeshMill.h"
#include "MeshMillDDG.h"
#include "MeshMapper.h"
#include "ShorelineMPI.h"
#include "MeshCourierServer.h"

// SCOREC includes
#include <apf.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <PCU.h>
#include <pcu_util.h>
#include <apfMDS.h>
#include <apfNumbering.h>
#include <apfNumberingClass.h>
#include <chef.h>
#include <phstream.h>
#include <phasta.h>

// vtk includes
#ifdef BUILD_WITH_VTK
  #include <vtkPolyData.h>
#endif

#ifdef USE_CATALYST
#include "vtkCPAdaptorAPI.h"
#include "vtkNew.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
//vtkCPProcessor* Processor = NULL;
void CatalystInit() {
/*
 * if(Processor == NULL) {
     Processor = vtkCPProcessor::New();
     Processor->Initialize();
   }
*/
  vtkCPAdaptorAPI::CoProcessorInitialize();
  vtkNew<vtkCPPythonScriptPipeline> pipeline;
  pipeline->Initialize("cpscript.py");
  vtkCPAdaptorAPI::GetCoProcessor()->AddPipeline(pipeline);
}
void CatalystFinalize() {
  vtkCPAdaptorAPI::CoProcessorFinalize();
  /*
   * if(Processor) {
   *   Processor->Delete();
   *   Processor = NULL;
   * }
  */
}
#endif

void setupChef(ph::Input& ctrl, int step) {
  //don't split or tetrahedronize
  ctrl.splitFactor = 1;
  ctrl.tetrahedronize = 0;
  ctrl.timeStepNumber = step;
  ctrl.solutionMigration = 1;
}

// colorize the output
const std::string reset("\033[0m");
const std::string red("\033[0;31m");
const std::string green("\033[0;32m");
const std::string magenta("\033[0;35m");

int main(int argc, char** argv) {
  // PetscInitialize calls MPI_Init
  PetscInitialize(&argc, &argv, (char*)0, NULL);
  PCU_Comm_Init();
  PCU_Protect();

#ifdef USE_CATALYST
  CatalystInit();
#endif

  // Shoreline::MeshCourier courier;
  // courier.CatalystInit();

  // Get the number of processes
  int worldSize;
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

  // Get the rank of the process
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  //////////////////////////////////////////////////////////////////////////////
  // read in mesh and geometry files using Chef
  gmi_model* g = nullptr;
  apf::Mesh2 * msh;
  gmi_register_mesh();
  gmi_register_null();

  if(!worldRank) std::cout << "setting up grs stream...\n";
  grstream grs = makeGRStream();

  // read in the Chef file adapt.inp
  ph::Input ctrl;
  if(!worldRank) std::cout << "loading adapt.inp...\n";
  ctrl.load("adapt.inp");

  // set some defaults
  ctrl.timeStepNumber = 0;   // set time step number to 0
  ctrl.writeRestartFiles = 0;// write a restart file to X_procs_case folder
  ctrl.writeGeomBCFiles = 1; // write a geombc file out to X_procs_case folder
  ctrl.solutionMigration = 0;// will crash if there is no IC (will set zero IC)
  ctrl.adaptFlag = 0;        // do not run mesh adaptation

  // load the mesh
  msh = apf::loadMdsMesh(ctrl.modelFileName.c_str(),ctrl.meshFileName.c_str());

  //////////////////////////////////////////////////////////////////////////////
  // setup numberings to be used by Shoreline
  msh->verify();
  //apf::reorderMdsMesh(msh);

  // setup mesh global numbering
  // apf::Sharing* shr = apf::getSharing(msh);
  apf::Numbering* ln = numberOwnedDimension(msh, "Local", 0);
  apf::Numbering* lnE = numberOwnedDimension(msh, "LocalElement", 3);
  apf::Numbering* lnOverlap = apf::numberOwnedNodes(msh, "Global");
  // apf::FieldShape* shape = getShape (ln);
  // apf::GlobalNumbering* gn = createGlobalNumbering (msh, "Global", shape);
  apf::GlobalNumbering* gn = apf::makeGlobal(lnOverlap);
  gn->rename("Global");
  synchronize(gn);

  msh->verify();

  //////////////////////////////////////////////////////////////////////////////
  // write geometry etc to the GRStream which will be passed to PHASTA
  chef::cook(g,msh,ctrl,grs);

  // define the RStream which serves as the output from PHASTA
  rstream rs = makeRStream();
  if(!worldRank) std::cout << "RStream created! reading PHASTA inputs\n";

  // read in the PHASTA input files solver.inp and input.config
  phSolver::Input inp("solver.inp", "input.config");

  //////////////////////////////////////////////////////////////////////////////
  // get a map to the boundary of the mesh
  Shoreline::MeshMapper mapper;

  //////////////////////////////////////////////////////////////////////////////
  // gather mesh onto a single process
  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();
  double before = MPI_Wtime();
  Shoreline::VectorMesh<long> vmParallel, vm0;
  std::vector<int> vertsPerProc;
  mapper.mapSurfaceToVectorMesh(msh, vmParallel);
  mapper.gatherVectorMesh(vmParallel, vm0, vertsPerProc, 0, MPI_COMM_WORLD);
  double after = MPI_Wtime();
  endTime = std::chrono::system_clock::now();
  if(!worldRank) {
    std::cout << "collected surface mesh on proc 0 in: " << (after - before)*1000 << " ms\n";
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //////////////////////////////////////////////////////////////////////////////
  // send VectorMesh to shClient
  std::vector<double> displacement;

  Shoreline::MeshServer server;
  if(worldRank == 0) {
    int port = 8080;
    server.attach(port);
    startTime = std::chrono::system_clock::now();
    bool didSendSurface = server.sendVectorMesh(vm0);
    endTime = std::chrono::system_clock::now();
    std::cout << "sent surface mesh to client in: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
              << " ms\n";
    if(didSendSurface)
      std::cout << "surface mesh sent with ";
    else
      std::cout <<"problem sending surface mesh with ";
    std::cout << vm0.vertIndices.size() << " vertices and "<< vm0.faceIndices.size()/3 << " faces\n";


#ifdef BUILD_WITH_VTK
    Shoreline::STLMesh<long> stlm;
    mapper.mapVectorMeshToSTLMesh(vm0, stlm);
    vtkSmartPointer<vtkPolyData> vtk = server.STLMesh2VTKPolyData(stlm);
    server.writeVTKPolyDataToFile(vtk, "serverSentMesh.vtp");
    std::cout << "server's surface mesh written to file\n";
#endif
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //////////////////////////////////////////////////////////////////////////////
  // time step!
  int step = 0;
  const int maxStep = 100000;
  const int deformationStart = 0;
  const int deformationStepCount = 1;
  const int deformationEnd = deformationStart + deformationStepCount;

  bool deformAvailable = false;
  int deformationSteps = 0;
  Shoreline::STLMesh<long> bcs;
  //////////////////////////////////////////////////////////////////////////////
  do {
    std::string vtkFilename = "./vtk/mesh" + std::to_string(step);
    if(step % 50 == 0) apf::writeASCIIVtkFiles(vtkFilename.c_str(), msh);

    if(!worldRank)
      std::cout << magenta << "phasta stepped to " << step << reset << "\n";

    // receive a displacement field
    std::cout << (!worldRank ? "checking for displacements...\n" : "");
    if(deformationSteps == 0 || deformationSteps >= deformationStepCount) { // check if we are done with last deformation
      bool receivedDisplacement = false;
      if(!worldRank) {
        receivedDisplacement = server.receiveVector(displacement, false);

      }
      // let all the other processes know we got a new displacement
      MPI_Bcast(&receivedDisplacement, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

      if(receivedDisplacement) {
        if(!worldRank) {
          if(displacement.size() != vm0.vertCoords.size()){
            std::cout << "displacement vector received was " << displacement.size() << " elements long\n";
            displacement.resize(vm0.vertCoords.size());
          }
          std::cout << "received displacement vector\n";
          std::cout << "displacement vector size: " <<displacement.size() << "\n";
          std::cout << "\n\n";
        }
        // scatter displacement field to all proceses
        const int fieldComponents = 3; //sending vector field
        std::vector<double> dispBCs;
        std::vector<int> sendCounts = vertsPerProc;
        for(auto& p : sendCounts){p *= fieldComponents;}
        Shoreline::scatterField(displacement, sendCounts, dispBCs, 0, MPI_COMM_WORLD);

        std::cout << "proc " << worldRank << " verts count = "
                  << vmParallel.vertCoords.size()/3 << ", nVerts = "
                  << vmParallel.nVerts << ", dispBC count = "
                  << dispBCs.size()/3 << "\n";
        vmParallel.vertCoords = dispBCs;
        std::cout << ((worldRank==0) ?"field scattered\n" : "");
        MPI_Barrier(MPI_COMM_WORLD);
        mapper.mapVectorMeshToSTLMesh(vmParallel, bcs);

        // linearly interpolate the deformation
        const double deformationScale = 1.0 / (deformationStepCount);
        for(auto & v : bcs.verts) {
          v.second[0] *= deformationScale;
          v.second[1] *= deformationScale;
          v.second[2] *= deformationScale;
        }

        // reset deformation step counter
        deformationSteps = 0;
	deformAvailable = true;
      } else {
        std::cout << ((worldRank==0) ? "no displacement field received\n": "");
      }
    }

    // deform mesh
    if( step >= deformationStart && deformAvailable) {
      deformationSteps++;
      if(deformationSteps >= deformationStepCount) deformAvailable = false;
      if(!worldRank)
        std::cout << green
                  << "starting mesh deformation at step "
                  << step << "..." << reset << "\n";
      before = MPI_Wtime();
      startTime = std::chrono::system_clock::now();
      Shoreline::MeshMill mill(msh);
      endTime = std::chrono::system_clock::now();
      after = MPI_Wtime();
      if(!worldRank) std::cout << "setup solver in " << (after - before)*1000 << " ms\n";

      before = MPI_Wtime();
      mill.formAndAssemble(bcs.verts);
      after = MPI_Wtime();
      if(!worldRank) std::cout << "form and assemble in " << (after - before)*1000 << " ms\n";

      before = MPI_Wtime();
      mill.solveMeshMotion();
      after = MPI_Wtime();
      if(!worldRank) std::cout << "form and assemble in " << (after - before)*1000 << " ms\n";
      mill.printPETScObjectInfo();
      mill.applyDeformation();
      mill.destroySolverObjects();
      endTime = std::chrono::system_clock::now();
      if(!worldRank){
        std::cout << "deformed mesh in: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
              << " ms\n";
        std::cout << green << "mesh deformed at step " << step << reset << "\n";
      }
    }

    step = phasta(inp,grs,rs);

    // apply the deformation after the flow solver runs
    clearGRStream(grs);
    setupChef(ctrl,step);
    chef::cook(g,msh,ctrl,rs,grs);
    clearRStream(rs);
  } while( step < maxStep );

  //////////////////////////////////////////////////////////////////////////////
  // cleanup
  MPI_Barrier(MPI_COMM_WORLD);
  // mill.destroySolverObjects();
  // destroy stream objects
  destroyGRStream(grs);
  destroyRStream(rs);
  // destroy mesh
  apf::destroyGlobalNumbering(gn);
  msh->destroyNative();
  apf::destroyMesh(msh);

  //////////////////////////////////////////////////////////////////////////////
  // Finalize the MPI environment.
  PCU_Comm_Free();

#ifdef USE_CATALYST
  CatalystFinalize();
#endif

  // turn off petsc and mpi
  PetscFinalize();
  return 0;
}
