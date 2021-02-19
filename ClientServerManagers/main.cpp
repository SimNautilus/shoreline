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

// vtk includes
#ifdef BUILD_WITH_VTK
  #include <vtkPolyData.h>
#endif
int main(int argc, char** argv){
  // PetscInitialize calls MPI_Init
  PetscInitialize(&argc, &argv, (char*)0, NULL);
  PCU_Comm_Init();
  // Shoreline::MeshCourier courier;
  // courier.CatalystInit();

  // Get the number of processes
  int worldSize;
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

  // Get the rank of the process
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  //////////////////////////////////////////////////////////////////////////////
  // read on the mesh. this will only work when running mpirun -np 4
  std::string geomDir = "../geometry/cube/";
  std::string modelFile = geomDir + "cube.dmg";
  std::string meshFile  = geomDir + "cube" + std::to_string(worldSize)
                        + "Par.smb";
  apf::Mesh2 * msh;
  gmi_register_mesh();
  gmi_register_null();
  msh = apf::loadMdsMesh(modelFile.c_str(), meshFile.c_str());
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

  //////////////////////////////////////////////////////////////////////////////
  // gather mesh onto a single process
  // int nVerts = 0, nFaces = 0, nVertsGlobal = 0, nFacesGlobal = 0;
  // std::map<long, std::array<double,3>> surfaceVerts;
  // std::vector<std::array<long,3>> surfaceConnectivity;
  // std::vector<long> vertIndsOnProc, faceIndsOnProc;
  // std::vector<double> vertCoordsOnProc;
  // mapper.mapSurfaceToVectorMesh(msh, nVerts, nFaces,
  //   vertIndsOnProc, vertCoordsOnProc, faceIndsOnProc);
  //
  // // get the total number of vertices and faces being sent to proc 0
  // MPI_Reduce(&nVerts, &nVertsGlobal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  // MPI_Reduce(&nFaces, &nFacesGlobal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  //
  // // gather all the vertices, vertex indices and face connectivity on proc 0
  // std::vector<long> vertIndices;
  // std::vector<double> vertCoords;
  // std::vector<long> faceIndices;
  // Shoreline::gatherVectors(vertIndsOnProc, MPI_LONG, vertIndices, 0 ,MPI_COMM_WORLD);
  // Shoreline::gatherVectors(vertCoordsOnProc, MPI_DOUBLE, vertCoords, 0 ,MPI_COMM_WORLD);
  // Shoreline::gatherVectors(faceIndsOnProc, MPI_LONG, faceIndices, 0 ,MPI_COMM_WORLD);
  // mapper.mapVectorMeshToSTLMesh(
  //   nVertsGlobal, nFacesGlobal,
  //   vertIndices, vertCoords, faceIndices,
  //   surfaceVerts, surfaceConnectivity
  // );

  //////////////////////////////////////////////////////////////////////////////
  // map received data on proc 0 to STL mesh and write to file
  // std::string filename = "../output/cube/srfMeshGlobal.vtp";;

  // for(int i = 0; i < 5; ++i) {
  //   courier.CatalystCoProcess(i, 0.0, surfaceVerts, surfaceConnectivity);
  //   cin.get();
  // }
  // if(worldRank == 0) {
  //   std::map<long, std::array<double,3>> surfaceVerts;
  //   std::vector<std::array<long,3>> surfaceConnectivity;
  //   // converting the data to a map removes duplicate entries from shared verts
  //   mapper.mapVectorMeshToSTLMesh(
  //     nVertsGlobal, nFacesGlobal,
  //     vertIndices, vertCoords, faceIndices,
  //     surfaceVerts, surfaceConnectivity
  //   );
  //   // map stl mesh to polydata mesh
  //   vtkSmartPointer<vtkPolyData> polydata =
  //     courier.STLMesh2VTKPolyData(surfaceVerts, surfaceConnectivity);
  //
  //   for(int i = 0; i < 5; ++i) {
  //     courier.CatalystCoProcess(i, 0.0, surfaceVerts, surfaceConnectivity);
  //     cin.get();
  //   }
  //   //write the polydata to a file
  //   courier.writeVTKPolyDataToFile(polydata, filename);
  //   std::cout << "global surface mesh written to file on process 0.\n";
  // }

  //////////////////////////////////////////////////////////////////////////////
  // create STLMesh on each process and print result to vtk file
  // std::map<long, std::array<double, 3>> verts;
  // std::vector<std::array<long, 3>> ien;
  // mapper.mapSurfaceToSTLMesh(msh, gn, verts, ien);

  // // map the STL mesh to a vtkPolyData object
  // vtkSmartPointer<vtkPolyData> polydata =
  //   courier.STLMesh2VTKPolyData(verts, ien);
  //
  // //write the polydata to a file
  // filename = "../output/cube/srfMesh_" + std::to_string(worldRank) + ".vtp";
  // courier.writeVTKPolyDataToFile(polydata, filename);
  //
  // std::cout << "Process " << worldRank
  //           << " wrote its piece of the surface mesh to a file.\n";

  //////////////////////////////////////////////////////////////////////////////
  // get boundary nodes on each process
  std::map<long, std::array<double, 3>> dispBCs;
  Shoreline::STLMesh<long> stlm;
  mapper.mapSurfaceToSTLMesh(msh, stlm);
  dispBCs = stlm.verts;

  //////////////////////////////////////////////////////////////////////////////
  // instantiate a Shoreline suite of objects
  Shoreline::MeshMill mill(msh);

  // just make up some boundary conditions
  for(auto & v : dispBCs) {
    double x = v.second[0];
    double y = v.second[1];
    double z = v.second[2];
    // double dispx = x * 0.25*exp(15 * y);
    // double dispy = .05 * y;
    // double dispz = 0.25 * z;

    double dispx = x * cos(10 * M_PI * y) - z * sin(10 * M_PI * y) - x;
    double dispy = 0;
    double dispz = x * sin(10 * M_PI * y) + z * cos(10 * M_PI * y) - z;
    v.second[0] = dispx;
    v.second[1] = dispy;
    v.second[2] = dispz;
  }

  // setup solver system
  mill.formAndAssemble(dispBCs);
  mill.solveMeshMotion();

  // apply deformation to the mesh
  mill.applyDeformation();
  // mill.viewSolution();
  // mill.viewMatrix();
  // write the deformed mesh to a file
  apf::writeASCIIVtkFiles("../output/cube/deformedMsh2", msh);

  //////////////////////////////////////////////////////////////////////////////
  // free the mesh
  MPI_Barrier(MPI_COMM_WORLD);
  apf::destroyGlobalNumbering(gn);
  // std::cout << "global numbering DESTROYED\n";
  mill.destroySolverObjects();
  // std::cout << "petsc objects DESTROYED\n";

  //////////////////////////////////////////////////////////////////////////////
  // Finalize the MPI environment.
  // courier.CatalystFinalize();
  PCU_Comm_Free();
  PetscFinalize(); // this should call MPI_Finalize
  // MPI_Finalize();
}
