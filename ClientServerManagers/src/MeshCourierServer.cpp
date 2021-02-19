
//stl includes
#include <iostream>
#include <fcntl.h>
#include <cassert>
// Shoreline includes
#include "MeshCourierServer.h"

namespace Shoreline {
bool MeshServer::attach(int port) {
  m_port = port;
  int opt = 1;
  size_t addrLen = sizeof(this->m_address);

  /* Creating sockets. AF_INET indicates that we are using an IPv4 domain
   * SOCK_STREAM "Provides sequenced, reliable, two-way, connection-
                  based byte streams.  An out-of-band data transmission
                  mechanism may be supported."
   * protocol is set to 0 because for basic stuff, the first protocol in the
   * domain (AF_INET) is all you need
   */
  if ((m_server = socket(AF_INET, SOCK_STREAM, 0)) == 0) {
      perror("socket failed");
      return false;
  }

  /* setsockopt provides a way of setting options on a socket
   * here, we are setting SO_REUSEPORT to true so we can reuse the port
   * SOL_SOCKET sets the option at the API level. If an option lives on
   * another level (like a TCP level) then you can specify which level to set
   * the option at.
   * SO_REUSEPORT allows multiple threads or processes on the same device
   * to also use the port. I suppose this will be useful for an hpc setting,
   * though for this toy problem, its not likely super important to set.
   */
  if (setsockopt(m_server, SOL_SOCKET, SO_REUSEPORT, &opt, sizeof(opt))) {
      perror("setsockopt");
      return false;
  }

  /* this section sets the server to non-blocking. i.e. if it doesn't hear
   * from a client application, it will continue runnning
   */
  //int flags = fcntl(m_server, F_GETFL, 0);
  //assert(flags != -1);
  //fcntl(m_server, F_SETFL, flags | O_NONBLOCK);
  
  // define the data address
  this->m_address.sin_family = AF_INET;
  this->m_address.sin_addr.s_addr = INADDR_ANY;
  this->m_address.sin_port = htons(port);

  /* In order to setup a server, we need to bind-listen-accept a socket
   * bind assigns an address to the socket
   * listen starts listening for connections to the server
   * accept accepts a connection to a client. this typically blocks the
   * servers progression until a client connects
   */
  if (bind(m_server, (struct sockaddr *)&this->m_address, addrLen) < 0) {
      perror("bind failed");
      return false;
  }
  std::cout << "bound\n";

  if (listen(m_server, 3) < 0) {
      perror("listen failed");
      return false;
  }
  std::cout << "listened\n";

  //this->m_socket = -1;
  this->m_socket = accept(m_server,(struct sockaddr*)&this->m_address,(socklen_t*)&addrLen);
  if (this->m_socket < 0){
    perror("accept failed");
    return false;
  }
  std::cout << "accepted\n";
  return true;
}



////////////////////////////////////////////////////////////////////////////////
template <typename T>
bool MeshServer::sendVectorMesh(const VectorMesh<T>& vm) {
  // ping the client to make sure it's ready to receive a vector mesh
  if(!this->sendPing(messageType::sendSrfMesh)) {
    std::cout << "problem with message type in receiveVectorMesh\n";
    return false;
  }

  // send the client how many vertices it should be expecting
  send(this->m_socket, &vm.nVerts, sizeof(int), 0);
  // make sure the client knows how many vertices there will be
  int nVertsTest = 0;
  recv(this->m_socket, &nVertsTest, sizeof(int), 0);
  if(nVertsTest != vm.nVerts) {
    std::cout << "problem with nVerts in sendVectorMesh\n";
    return false;
  }

  // send the vertex coordinates
  size_t sendBytes = 0;
  while(sendBytes < 3 * vm.nVerts * sizeof(double)) {
    sendBytes += send(this->m_socket, (char*)vm.vertCoords.data() + sendBytes,
      3 * vm.nVerts * sizeof(double) - sendBytes, 0);
  }


  // send the vertex indices data
  sendBytes = 0;
  while(sendBytes < vm.nVerts * sizeof(T)) {
    sendBytes += send(this->m_socket, (char*)vm.vertIndices.data() + sendBytes,
    vm.nVerts * sizeof(T) - sendBytes, 0);
  }

  // send the client how many vertices it should be expecting
  send(this->m_socket, &vm.nFaces, sizeof(int), 0);
  // make sure the client knows how many vertices there will be
  int nFacesTest = 0;
  recv(this->m_socket, &nFacesTest, sizeof(int), 0);
  if(nFacesTest != vm.nFaces) {
    std::cout << "problem with nFaces in sendVectorMesh\n";
    return false;
  }

  // send the face indices data
  sendBytes = 0;
  while(sendBytes < 3 * vm.nFaces * sizeof(T)) {
    sendBytes += send(this->m_socket, (char*)vm.faceIndices.data() + sendBytes,
      3 * vm.nFaces * sizeof(T) - sendBytes, 0);
  }

  // send the featureTags data
  sendBytes = 0;
  while(sendBytes < vm.featureTags.size() * sizeof(int)) {
    sendBytes += send(this->m_socket, (char*)vm.featureTags.data() + sendBytes,
    vm.featureTags.size() * sizeof(int) - sendBytes, 0);
  }

  // exit hoping everything went ok
  return true;
}

// forward declarations
template bool MeshServer::sendVectorMesh<int>(const VectorMesh<int>& vm);
template bool MeshServer::sendVectorMesh<long>(const VectorMesh<long>& vm);

}
