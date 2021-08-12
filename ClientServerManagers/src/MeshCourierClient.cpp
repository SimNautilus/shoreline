#include "MeshCourierClient.h"
#include <chrono>

namespace Shoreline {
bool MeshClient::attach(const int port) {
  this->m_port = port;
  // create a couple sockets where we can connect to a server
  // kill the client if we cannot create them
  if ((this->m_socket = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
    std::cout << "\n Socket creation error \n";
    return false;
  }
  /* set the protocol for the network. AF_INET sets it to ipv4 connections
   * this could also be AF_INET6 for ipv7 or AF_BLUETOOTH for bluetooth
   * connections if you were into that.
  */
  this->m_address.sin_family = AF_INET;
  // converts a uint16_t from host byte ordering to network byte ordering
  this->m_address.sin_port = htons(port);
  // get the size of the address struct which we will need when we connect
  // to the socket
  size_t addrLen = sizeof(this->m_address);

  // Convert IPv4 and IPv6 addresses from text to binary form
  // the ip address 127.0.0.1 is special and indicates a localhost
  // the assumption here is that the user will connect to the server over
  // ssh -L<clientPort>:<serverhost>:<serverPort>
  // so the localhost will connect to the clientPort. hopefully we don't need
  // to change this ip address situation here
  if(inet_pton(AF_INET, "127.0.0.1", &this->m_address.sin_addr)<=0) {
    std::cout << "\nInvalid address / Address not supported \n";
    return false;
  }

  // actually connect to the socket
  if (connect(this->m_socket,(struct sockaddr *)&this->m_address, addrLen) < 0) {
    std::cout << "\nConnection Failed \n";
    return false;
  } else {
    std::cout << "\nConnection Established!\n";
    return true;
  }
}

////////////////////////////////////////////////////////////////////////////////
template <typename T>
bool MeshClient::receiveVectorMesh(VectorMesh<T>& vm) {
  // wait for ping from server that its going to send a vector mesh
  // this is a blocking communication (hopefully)
  if(!this->receivePing(messageType::sendSrfMesh)) {
    std::cout << "problem with message type in receiveVectorMesh\n";
    return false;
  }

  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();
  
  // wait for info on how many vertices we will be receiving
  recv(this->m_socket, &vm.nVerts, sizeof(int), 0);
  std::cout << "nVerts = " << vm.nVerts <<"\n";
  // let the server know we got the verts
  send(this->m_socket, &vm.nVerts, sizeof(int), 0);

  // resize vm.vertCoords and vm.vertIndices
  vm.vertCoords.resize(3 * vm.nVerts);
  vm.vertIndices.resize(vm.nVerts);

  // receive the vertex coordinates
  size_t recvBytes = 0;
  while(recvBytes < 3 * vm.nVerts * sizeof(double)) {
    recvBytes += recv(this->m_socket, (char*)vm.vertCoords.data()+recvBytes,
      3 * vm.nVerts * sizeof(double) - recvBytes, 0);
  }


  // receive the vertex indices data
  recvBytes = 0;
  while(recvBytes < vm.nVerts * sizeof(T)) {
    recvBytes += recv(this->m_socket, (char*)vm.vertIndices.data() + recvBytes,
      vm.nVerts * sizeof(T) - recvBytes, 0);
  }



  // wait for info on how many faces we will be receiving
  recv(this->m_socket, &vm.nFaces, sizeof(int), 0);
  // let the server know we got the faces
  send(this->m_socket, &vm.nFaces, sizeof(int), 0);

  // resize the vm.faceIndices
  vm.faceIndices.resize(3 * vm.nFaces);
  // receive the face indices data
  recvBytes = 0;
  while(recvBytes < 3 * vm.nFaces * sizeof(T)) {
    recvBytes += recv(this->m_socket, (char*)vm.faceIndices.data()+recvBytes,
      3 * vm.nFaces * sizeof(T) - recvBytes, 0);
  }

  // receive the featureTags data
  vm.featureTags.resize(vm.nFaces);
  recvBytes = 0;
  while(recvBytes < vm.nFaces * sizeof(int)) {
    recvBytes += recv(this->m_socket, (char*)vm.featureTags.data() + recvBytes,
      vm.nFaces * sizeof(int) - recvBytes, 0);
  }

  endTime = std::chrono::system_clock::now();
  std::cout << "received vector mesh in :"
            <<std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
            << " ms\n";
  // exit hoping everything went ok
  return true;
}

// forward declarations
template bool MeshClient::receiveVectorMesh<int>(VectorMesh<int>& vm);
template bool MeshClient::receiveVectorMesh<long>(VectorMesh<long>& vm);

}
