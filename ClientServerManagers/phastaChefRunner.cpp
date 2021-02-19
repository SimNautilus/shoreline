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
#include "MeshCourier.h"
#include "MeshMapper.h"
#include "ShorelineMPI.h"

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

void setupChef(ph::Input& ctrl, int step) {
  //don't split or tetrahedronize
  ctrl.splitFactor = 1;
  ctrl.tetrahedronize = 0;
  ctrl.timeStepNumber = step;
  ctrl.solutionMigration = 1;
}


int main(int argc, char** argv){
  // PetscInitialize calls MPI_Init
  PetscInitialize(&argc, &argv, (char*)0, NULL);
  PCU_Comm_Init();
  PCU_Protect();
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
  apf::reorderMdsMesh(msh);
  
  // setup mesh global numbering
  apf::Sharing* shr = apf::getSharing(msh);
  apf::Numbering* ln = numberOwnedDimension(msh, "Local", 0);
  apf::Numbering* lnE = numberOwnedDimension(msh, "LocalElement", 3);
  apf::Numbering* lnOverlap = apf::numberOwnedNodes(msh, "Global", 0);

  apf::GlobalNumbering* gn = apf::makeGlobal(lnOverlap);
  gn->rename("Global");
  synchronize(gn, shr);

  //////////////////////////////////////////////////////////////////////////////
  // write geometry etc to the GRStream which will be passed to PHASTA
  chef::cook(g,msh,ctrl,grs);
  if(!worldRank) std::cout << "chef has been cooked! defining RStream...\n";
  char** ggptr = (char**)grs;
  char* maybeRestart = ggptr[0];
  std::cout << maybeRestart << "\n";
  std::cout << sizeof(size_t) << "\n";
  std::cout << "size of gstream gSz = " << (size_t)ggptr[2] << " size of restart rSz = " << (size_t)ggptr[3] << "\n";
  if(ggptr[0] == NULL){std::cout << "\nWE GOT A NULLPTR ON THE geombc PROBABLY\n";}
  if(ggptr[1] == NULL){std::cout << "\nWE GOT A NULLPTR ON THE RESTART PROBABLY\n";}
  // define the RStream which serves as the output from PHASTA
  rstream rs = makeRStream();
  if(!worldRank) std::cout << "RStream created! reading PHASTA inputs\n";

  // read in the PHASTA input files solver.inp and input.config
  phSolver::Input inp("solver.inp", "input.config");
  
  //////////////////////////////////////////////////////////////////////////////
  // get a map to the boundary of the mesh
  Shoreline::MeshMapper mapper;

  // get boundary nodes on each process
  Shoreline::STLMesh<long> stlm;
  mapper.mapSurfaceToSTLMesh(msh, stlm);

  //////////////////////////////////////////////////////////////////////////////
  // instantiate a Shoreline suite of objects
  // apf::writeASCIIVtkFiles("../output/cubeBL/undeformedMsh", msh);
  Shoreline::MeshMill mill(msh);

  // just make up some boundary conditions
  for(auto & v : stlm.verts) {
    double x = v.second[0];
    double y = v.second[1];
    double z = v.second[2];

    // rotate about x axis
    double dispx = 0;
    double dispy = (y-0.05) * cos(5 * M_PI * (x+0.05)) 
                 - z * sin(5 * M_PI * (x+0.05)) - y + 0.05;
    double dispz = (y-0.05) * sin(5 * M_PI * (x+0.05)) 
                 + z * cos(5 * M_PI * (x+0.05)) - z;
    v.second[0] = dispx;
    v.second[1] = dispy;
    v.second[2] = dispz;
  }

  // setup solver system
  using namespace std::chrono;
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  mill.formAndAssemble(stlm.verts);
  mill.printPETScObjectInfo();
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

  if(!worldRank)std::cout<<"Form and Assemble: "<<time_span.count()<<" seconds.\n";

  // solve for the mesh deformation (do not apply yet)
  t1 = high_resolution_clock::now();
  mill.solveMeshMotion();
  mill.printPETScObjectInfo();
  t2 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t2 - t1);
  if(!worldRank) std::cout << "Solve: " << time_span.count() << " seconds.\n";

  //////////////////////////////////////////////////////////////////////////////
  // time step!
  int step = 0;
  int maxStep = 2;
  do {
    std::string vtkFilename = "./vtk/mesh" + std::to_string(step);

    if(!worldRank) std::cout << "about to phasta step...\n";
    step = phasta(inp,grs,rs);
    if(!worldRank) std::cout << "phasta Stepped!\n";
    
    clearGRStream(grs);
    if(!PCU_Comm_Self()) fprintf(stderr, "CAKE ran to step %d\n", step);
    setupChef(ctrl,step);
//    mill.applyDeformation();
    chef::cook(g,msh,ctrl,rs,grs);
    clearRStream(rs);
    apf::writeASCIIVtkFiles(vtkFilename.c_str(), msh);
    // apply the deformation after the flow solver runs
  } while( step < maxStep );
 
  //////////////////////////////////////////////////////////////////////////////
  // cleanup
  MPI_Barrier(MPI_COMM_WORLD);

  // destroy the 
  apf::destroyGlobalNumbering(gn);
  mill.destroySolverObjects();
  destroyGRStream(grs);
  destroyRStream(rs);
  msh->destroyNative();
  apf::destroyMesh(msh);

  //////////////////////////////////////////////////////////////////////////////
  // Finalize the MPI environment.
  // courier.CatalystFinalize();
  PCU_Comm_Free();
  PetscFinalize(); // this should call MPI_Finalize
}
