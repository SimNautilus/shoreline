/* @file MeshClient.h
 * @brief Header file for the MeshClient class
 * @details MeshClient provides routines for receiving data from a MeshServer
 *          object.
 *
 */

 #pragma once

// Shoreline includes
#include "MeshCourier.h"
#include "MeshMapper.h"

namespace Shoreline {
class MeshClient : public MeshCourier {
public:
  // @brief default constructor
  MeshClient() : MeshCourier(){};

  // @brief attach to a server with the given port number
  bool attach(const int port) override;

  // @brief receive a VectorMesh from a connected server
  template <typename T>
  bool receiveVectorMesh(VectorMesh<T>& vm);

protected:
};
}
