#include "MeshMapper.h"
#include "MeshCourier.h"
#include "MeshCourierClient.h"

#include <iostream>
#include <fstream>
#include <chrono>

int main() {
  Shoreline::VectorMesh<long> vm;
  Shoreline::MeshClient client;

  // receive VectorMesh from server
  int port = 8080;
  client.attach(port);
  client.receiveVectorMesh(vm);
  std::cout << "mesh received. NVerts = "
    << vm.nVerts << " nFaces = " << vm.nFaces << "\n";

  //////////////////////////////////////////////////////////////////////////////
  Shoreline::MeshMapper mapper;
  Shoreline::STLMesh<long> stlm;
  mapper.mapVectorMeshToSTLMesh(vm, stlm);

  //////////////////////////////////////////////////////////////////////////////
  // write VTK file of the surface mesh
  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();
  #ifdef BUILD_WITH_VTK
    vtkSmartPointer<vtkPolyData> vtk = client.STLMesh2VTKPolyData(stlm);
    client.writeVTKPolyDataToFile(vtk, "clientReceivedMesh.vtp");
    std::cout << "mesh boundary written to vtp file\n";
  #endif
  endTime = std::chrono::system_clock::now();
  std::cout << "mesh exported to ParaView in :"
            <<std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
            << " ms\n";
  //////////////////////////////////////////////////////////////////////////////
  // hang out and wait for ParaView to deal with modifying the mesh
  std::vector<double> displacement;
  displacement.reserve(3 * stlm.verts.size());

  bool quit = false;
  while(!quit) {
    // wait for the user to hit enter
    std::cout << "waiting for user's signal that displacment.csv is "
              << "written to current working directory...\n"
              << "  press \'s\' to confirm displacement field\n"
              << "  press \'q\' to quit\n";

    char input = std::cin.get();

    if(input == 's' || input == 'S') {
      //when the user has notified, check the filepath for a displacement.txt
      std::fstream dispFile("displacements.csv", std::ios_base::in);
      if(!dispFile.is_open()) {
        std::cout << "PROBLEM OPENING FILE\n";
      }
      double xx = 0;
      displacement.resize(0);
      while(dispFile >> xx) {
        displacement.push_back(xx);
      }
      std::cout << "displacement vector size: " <<displacement.size() << "\n";
      if(displacement.size() < 3 * stlm.verts.size()) {
        displacement.resize(3 * stlm.verts.size(), 0.0);
      }

      ////////////////////////////////////////////////////////////////////////////
      // inflate the displacement field back to a distributable size (with overlaps in place)
      std::vector<double> displacementInflated(vm.vertCoords.size());
      for(int i = 0; i < vm.nVerts; ++i) {
        auto it = stlm.verts.find(vm.vertIndices[i]);
        int idisp = std::distance(stlm.verts.begin(),it);
        displacementInflated[3*i  ] = displacement[3*idisp  ];
        displacementInflated[3*i+1] = displacement[3*idisp+1];
        displacementInflated[3*i+2] = displacement[3*idisp+2];
      }
      std::cout << "inflated displacement vector size: " <<displacementInflated.size() << "\n";
      // send the displacement field back to the server
      client.sendVector(displacementInflated);
    }
    if(input == 'q' || input == 'Q') {
      quit = true;
    }
  }

  return 0;
}
