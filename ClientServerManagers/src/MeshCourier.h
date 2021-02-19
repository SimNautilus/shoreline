/* @file MeshCourier.h
 * @brief Header file for the MeshCourier class
 * @details MeshCourier provides tools for passing a surface mesh to and from
 *          a user for modification on a user device. This is done by taking a
 *          map to a surface mesh (stored distributed across all processes) and
 *          writing a vtkPolyData object which is then sent to the user over a
 *          network interconnect managed by <CATALYST>
 *          (probably. Probably not Sensei... TBD). The user can then modify
 *          the mesh on their local computer and send the modified surface mesh
 *          back to the MeshCourier object residing on the compute resouce.
 *          This modified surface mesh is then converted distributed to all
 *          cores on the compute resource and that data can be used by the
 *          MeshMill as boundary conditions to compute a volumetric deformation
 *          of the mesh domain.
 *
 */
#pragma once
// vtk includes
#ifdef BUILD_WITH_VTK
  #include <vtkSmartPointer.h>
  #include <vtkPolyData.h>
  #include <vtkTriangle.h>
  #include <vtkCellArray.h>
#endif
// stl includes
#include <map>
#include <vector>
#include <array>

// Shoreline Includes
#include "MeshMapper.h"

// network system includes
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
////////////////////////////////////////////////////////////////////////////////
namespace Shoreline {
class MeshCourier {
public:
  MeshCourier() {
    m_port = 0;
    m_socket = 0;
  };
  void print();

  enum messageType {
    nullMsg = -1,
    sendVecMsg,
    sendSrfMesh,
    sendDeformationField
  };

  // pure virtual function attach to another MeshCourier derived object
  virtual bool attach(const int port) = 0;

  // @brief ping the connected courier with a message type
  bool sendPing(messageType msg);

  // @brief receive ping from connected courier
  bool receivePing(messageType msg);

  // @brief send a vector to the server (could be displacement fields etc.)
  template <typename T> bool sendVector(const std::vector<T>& vec);

  // @brief receive vector from client (displacement fields etc.)
  template <typename T> bool receiveVector(std::vector<T>& vec, bool waitForData);

  /* @brief write a vtkPolyData object to a file given the filename */
#ifdef BUILD_WITH_VTK
  static void writeVTKPolyDataToFile(
    vtkSmartPointer<vtkPolyData> polydata, std::string filename
  );

  /* @brief convert an STL mesh to a vtkPolyData */
  template <typename T>
  static vtkSmartPointer<vtkPolyData> STLMesh2VTKPolyData(STLMesh<T>&stlm);
#endif

protected:
  // port number to connect to
  int m_port = 0;

  // @brief socket descriptor
  int m_socket = 0;

  // @brief struct to store info about the socket
  struct sockaddr_in m_address;

private:
};
} // end namespace Shoreline
