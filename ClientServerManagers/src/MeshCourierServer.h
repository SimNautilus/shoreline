/* @file MeshServer.h
 * @brief Header file for the MeshServer class
 * @details MeshServer provides routines for sending data to a MeshClient
 *          object.
 *
 */

#pragma once

// Shoreline includes
#include "MeshCourier.h"
#include "MeshMapper.h"


namespace Shoreline {
class MeshServer : public MeshCourier {
public:
  // @brief default constructor
  MeshServer() : MeshCourier(){};

  // @brief attach to a client with the given port number
  bool attach(int port) override;

  // @brief send a VectorMesh to a connected client
  template <typename T>
  bool sendVectorMesh(const VectorMesh<T>& vm);

protected:
  // @brief server descriptor
  int m_server = 0;

};
}
