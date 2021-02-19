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



// colorize the output
const std::string reset("\033[0m");
const std::string red("\033[0;31m");
const std::string green("\033[0;32m");
const std::string magenta("\033[0;35m");
const std::string yellow("\033[0;33m");

int main(int argc, char** argv) {
  // PetscInitialize calls MPI_Init
  PetscInitialize(&argc, &argv, (char*)0, NULL);
  PCU_Comm_Init();
  PCU_Protect();

  std::string vtkFilename = "msh";
  int maxStep = 1;
  if(argc > 1)
    vtkFilename = argv[1];
  if(argc > 2)
    maxStep = atoi(argv[2]);

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

  // load the mesh
  msh = apf::loadMdsMesh(ctrl.modelFileName.c_str(),ctrl.meshFileName.c_str());

  //////////////////////////////////////////////////////////////////////////////
  // setup numberings to be used by Shoreline
  msh->verify();
  apf::reorderMdsMesh(msh);

  // set some defaults
  ctrl.splitFactor = 1;
  ctrl.tetrahedronize = 0;
  ctrl.timeStepNumber = 0;
  ctrl.writeRestartFiles = 0;
  ctrl.writeGeomBCFiles = 1;
  ctrl.solutionMigration = 0;
  ctrl.adaptFlag = 0;

  if(!worldRank) std::cout << "adapt.inp loaded! cooking chef...\n";
  chef::preprocess(msh,ctrl,grs);
  // chef::cook(g,msh,ctrl,grs);
  if(!worldRank) std::cout << "chef has been cooked! defining RStream...\n";
  rstream rs = makeRStream();
  if(!worldRank) std::cout << "RStream created!\n";

  // setup mesh global numbering
  apf::Numbering* ln = numberOwnedDimension(msh, "Local", 0);
  apf::Numbering* lnE = numberOwnedDimension(msh, "LocalElement", 3);
  apf::Numbering* lnOverlap = apf::numberOwnedNodes(msh, "Global");

  apf::GlobalNumbering* gn = apf::makeGlobal(lnOverlap);
  gn->rename("Global");
  synchronize(gn);

  msh->verify();
  //////////////////////////////////////////////////////////////////////////////
  // get a map to the boundary of the mesh
  Shoreline::MeshMapper mapper;

  //////////////////////////////////////////////////////////////////////////////
  // gather mesh onto a single process
  Shoreline::VectorMesh<long> vmParallel, vm0;
  std::vector<int> vertsPerProc;
  mapper.mapSurfaceToVectorMesh(msh, vmParallel);
  mapper.gatherVectorMesh(vmParallel, vm0, vertsPerProc, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  //////////////////////////////////////////////////////////////////////////////
  // send VectorMesh to shClient
  std::vector<double> displacement;

  if(worldRank == 0) {
    int port = 8080;
    Shoreline::MeshServer server;
    server.attach(port);
    server.sendVectorMesh(vm0);
    std::cout << "surface mesh sent\n";

    // receive a displacement field
    std::cout << "waiting for displacement vector...\n";
    if(server.receiveVector(displacement, true))
      std::cout << "received displacement vector\n";
    // apply displacement field
    std::cout << "displacement vector size: " <<displacement.size() << "\n";
    std::cout << "\n\n";
    if(displacement.size() != vm0.vertCoords.size()) {
      std::cout << "resizing displacement vector to "
                << vm0.vertCoords.size() << "\n";
      displacement.resize(vm0.vertCoords.size(),0.0);
    }
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
  Shoreline::STLMesh<long> bcs;
  std::cout << ((worldRank==0) ?"field scattered\n" : "");
  MPI_Barrier(MPI_COMM_WORLD);
  mapper.mapVectorMeshToSTLMesh(vmParallel, bcs);

  //////////////////////////////////////////////////////////////////////////////
  // perform mesh modification on the volume mesh
  int step = 0;
  Shoreline::STLMesh<long> bcScaled = bcs;

  // linearly interpolate the deformation
  const double deformationScale = 1.0/maxStep;
  for(auto & v : bcs.verts){
    v.second[0] *= deformationScale;
    v.second[1] *= deformationScale;
    v.second[2] *= deformationScale;
  }

  do {
    if(!worldRank)
      std::cout << yellow << "starting step " << step << "/" << maxStep << reset << "\n";
    Shoreline::MeshMill mill(msh);
    mill.formAndAssemble(bcs.verts);
    mill.solveMeshMotion();
    mill.printPETScObjectInfo();
    mill.applyDeformation();
    msh->verify();

    apf::writeASCIIVtkFiles(vtkFilename.c_str(), msh);
    mill.destroySolverObjects();
    if(!worldRank)
      std::cout << green << "mesh deformed at step " << step << "/" << maxStep << reset << "\n";
    step++;
  } while( step < maxStep );

  //////////////////////////////////////////////////////////////////////////////
  // cleanup
  MPI_Barrier(MPI_COMM_WORLD);

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
  // turn off petsc and mpi
  PetscFinalize();
  return 0;
}
