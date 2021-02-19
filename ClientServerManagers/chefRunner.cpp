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


// vtk includes
#ifdef BUILD_WITH_VTK
  #include <vtkPolyData.h>
#endif

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
  gmi_model* g = 0;
  apf::Mesh2 * msh;
  gmi_register_mesh();
  gmi_register_null();

  std::cout << "setting up grs stream...\n";
  grstream grs = makeGRStream();
  ph::Input ctrl;
  std::cout << "loading adapt.inp...\n";
  ctrl.load("adapt.inp");
  msh = apf::loadMdsMesh(
      ctrl.modelFileName.c_str(),ctrl.meshFileName.c_str()
  );
  // set some defaults
  ctrl.splitFactor = 1;
  ctrl.tetrahedronize = 0;
  ctrl.timeStepNumber = 0;
  ctrl.writeRestartFiles = 1;
  ctrl.writeGeomBCFiles = 1;
  ctrl.solutionMigration = 0;
  ctrl.adaptFlag = 0;

  std::cout << ctrl.meshFileName << ", " << ctrl.modelFileName << "\n";
  std::cout << "adapt.inp loaded! cooking chef...\n";
  chef::preprocess(msh,ctrl,grs);
  // chef::cook(g,msh,ctrl,grs);
  std::cout << "chef has been cooked! defining RStream...\n";
  rstream rs = makeRStream();
  std::cout << "RStream created!\n";

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
  // get a map to the boundary of the mesh
  Shoreline::MeshMapper mapper;

  // get boundary nodes on each process
  std::map<long, std::array<double, 3>> dispBCs;
  Shoreline::STLMesh<long> stlm;
  mapper.mapSurfaceToSTLMesh(msh, stlm);
  dispBCs = stlm.verts;

  //////////////////////////////////////////////////////////////////////////////
  // instantiate a Shoreline suite of objects
  // apf::writeASCIIVtkFiles("../output/cubeBL/undeformedMsh", msh);
  Shoreline::MeshMill mill(msh);

  // just make up some boundary conditions
  for(auto & v : dispBCs) {
    double x = v.second[0];
    double y = v.second[1];
    double z = v.second[2];

    // rotate about x axis
    double dispx = 0;
    double dispy = (y-0.05) * cos(5 * M_PI * (x+0.05)) - z * sin(5 * M_PI * (x+0.05)) - y + 0.05;
    double dispz = (y-0.05) * sin(5 * M_PI * (x+0.05)) + z * cos(5 * M_PI * (x+0.05)) - z;
    v.second[0] = dispx;
    v.second[1] = dispy;
    v.second[2] = dispz;
  }

  // setup solver system
  using namespace std::chrono;
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  mill.formAndAssemble(dispBCs);
  mill.printPETScObjectInfo();
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

  std::cout << "Form and Assemble: " << time_span.count() << " seconds.\n";

  t1 = high_resolution_clock::now();
  mill.solveMeshMotion();
  mill.printPETScObjectInfo();
  t2 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t2 - t1);
  std::cout << "Solve: " << time_span.count() << " seconds.\n";

  // apply deformation to the mesh
  mill.applyDeformation();
  // mill.viewSolution();
  // mill.viewMatrix();
  // write the deformed mesh to a file
  apf::writeASCIIVtkFiles("output/cubeBL/deformedMshLE_Serial", msh);

  //////////////////////////////////////////////////////////////////////////////
  // free the mesh
  MPI_Barrier(MPI_COMM_WORLD);
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
  // MPI_Finalize();
}
