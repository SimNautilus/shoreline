#include "MeshCourier.h"
#include <iostream>

// vtk includes
#ifdef BUILD_WITH_VTK
  #include <vtkXMLPolyDataWriter.h>
  #include <vtkDoubleArray.h>
  #include <vtkCellData.h>
#endif


namespace Shoreline {
////////////////////////////////////////////////////////////////////////////////
bool MeshCourier::sendPing(messageType msg) {
  send(m_socket, &msg, sizeof(int), 0);
  // wait for response ping to make sure we are on the same page as the client
  messageType msgRecv = nullMsg;
  recv(m_socket, &msgRecv, sizeof(int), 0);
  if(msg != msgRecv) {
    return false;
  }
  return true;
}

bool MeshCourier::receivePing(messageType msg) {
  messageType msgRecv = messageType::nullMsg;
  //recv(m_socket, &msgRecv, sizeof(int), MSG_PEEK | MSG_DONTWAIT);
  //if(msgRecv != msg) return false;

  recv(m_socket, &msgRecv, sizeof(int), 0);
  // let the client know we are on the same page
  send(m_socket, &msg, sizeof(int), 0);
  if(msg != msgRecv) {
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
template <typename T>
bool MeshCourier::receiveVector(std::vector<T>& vec, bool waitForData) {
  if(waitForData) {
    while(!receivePing(messageType::sendVecMsg)){}
  } else {
    if(!receivePing(messageType::sendVecMsg)) {return false;}
  }

  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  // receive the number of elements in the vector
  int nElements = 0;
  recv(this->m_socket, &nElements, sizeof(int), 0);
  // parrot the number of elements sent
  send(this->m_socket, &nElements, sizeof(int), 0);

  // resize vector to receive data
  vec.resize(nElements);

  startTime = std::chrono::system_clock::now();
  // receive the data
  size_t recvBytes = 0;
  while(recvBytes < vec.size() * sizeof(T)) {
    recvBytes += recv(this->m_socket, (char*)vec.data()+recvBytes,
      vec.size() * sizeof(T) - recvBytes, 0);
  }
  endTime = std::chrono::system_clock::now();
  std::cout << "received vector in :" <<std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms\n";
  return true;
}

// template specializations for several useful types
template bool MeshCourier::receiveVector<int>(std::vector<int>& vec, bool waitForData = true);
template bool MeshCourier::receiveVector<long>(std::vector<long>& vec, bool waitForData = true);
template bool MeshCourier::receiveVector<double>(std::vector<double>& vec, bool waitForData = true);

////////////////////////////////////////////////////////////////////////////////
template <typename T>
bool MeshCourier::sendVector(const std::vector<T>& vec) {
  while(!sendPing(messageType::sendVecMsg)) {}

  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;

  startTime = std::chrono::system_clock::now();

  // send the number of elements in the vector
  int nElements = vec.size();
  send(this->m_socket, &nElements, sizeof(int), 0);

  // check that the other courier is ready to receive data
  int nElmsTest = 0;
  recv(this->m_socket, &nElmsTest, sizeof(int), 0);
  if(nElmsTest != nElements) {
    std::cout << "problem with nElements in sendVector\n";
    return false;
  }

  // send the vector data
  size_t sentBytes = 0;
  while(sentBytes < vec.size() * sizeof(T)) {
    sentBytes += send(this->m_socket, (char*)vec.data() + sentBytes,
      vec.size() * sizeof(T)-sentBytes, 0);
  }


  endTime = std::chrono::system_clock::now();
  std::cout << "sent vector in :" <<std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms";
  return true;
}

// template specializations for several useful types
template bool MeshCourier::sendVector<int>(const std::vector<int>& vec);
template bool MeshCourier::sendVector<long>(const std::vector<long>& vec);
template bool MeshCourier::sendVector<double>(const std::vector<double>& vec);

////////////////////////////////////////////////////////////////////////////////
void MeshCourier::print(){
  std::cout << "MeshCourier Instanced\n";
}

////////////////////////////////////////////////////////////////////////////////
#ifdef BUILD_WITH_VTK
void MeshCourier::writeVTKPolyDataToFile(
  vtkSmartPointer<vtkPolyData> polydata, std::string filename)
{
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName (filename.c_str());
  writer->SetInputData ( polydata );
  writer->Write();
}

template <typename T>
vtkSmartPointer<vtkPolyData> MeshCourier::STLMesh2VTKPolyData(STLMesh<T>&stlm)
{
  // extract the points into a vtk array
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  for(auto & ptMap : stlm.verts) {
    auto pt = ptMap.second;
    points->InsertNextPoint(pt[0], pt[1], pt[2]);
  }

  // extract the triangles
  vtkSmartPointer<vtkCellArray> tris = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
  vtkSmartPointer<vtkDoubleArray> tags = vtkSmartPointer<vtkDoubleArray>::New();
  tags->SetName("features");
  tags->SetNumberOfComponents(1);
  tags->SetNumberOfTuples(stlm.featureTags.size());

  for(auto & tri : stlm.ien) {
    for(int i = 0; i < 3; ++i) {
      // find the location of the vertex in the map O(log n) lookup
      auto it = stlm.verts.find(tri[i]);
      unsigned int index = std::distance(stlm.verts.begin(), it);
      triangle->GetPointIds()->SetId(i, index);
    }
    tris->InsertNextCell(triangle);
  }

  // wite a feature tag list to the vtkpolydata
  unsigned int index = 0;
  for(auto tag : stlm.featureTags) {
    tags->SetValue(index++, tag);
  }

  // define a vtkPolyData
  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();

  polydata->SetPoints(points);
  polydata->SetPolys(tris);
  polydata->GetCellData()->AddArray(tags);
  return polydata;
}

template vtkSmartPointer<vtkPolyData> MeshCourier::STLMesh2VTKPolyData<int>(STLMesh<int>&stlm);
template vtkSmartPointer<vtkPolyData> MeshCourier::STLMesh2VTKPolyData<long>(STLMesh<long>&stlm);
#endif
} // end namespace Shoreline
