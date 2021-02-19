/* @file MeshMapper.h
 * @brief Header file for the MeshMapper class
 * @details MeshMapper provides tools for generating maps across aspects of a
 *          mesh domain. This includes:
 *          - A surface mesh map which can be passed to a MeshCourier
 *            and carted off to a user for modification
 *
 */
#pragma once

// std includes
#include <vector>
#include <array>
#include <map>
#include <iostream>
#include <mpi.h>
#include <algorithm>

// SCOREC includes
#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <apfNumbering.h>
#include <apfShape.h>

namespace Shoreline {

/* @brief VectorMesh provides a simple data structure for passing around a mesh
 * @note currently only supports triangulated surfaces
 * @note templated on index integer size (int, long int, uint, etc.)
*/
template<typename T>
struct VectorMesh {
  int nVerts;
  int nFaces;
  std::vector<T> vertIndices; // size nVerts
  std::vector<double> vertCoords; // size 3*nVerts [x1,y1,z1,...,xn,yn,zn]
  std::vector<T> faceIndices; // size 3*Faces [v11,v12,v13,...,vn1,vn2,vn3]
  std::vector<int> featureTags; // geometric tags for each vertex
};

/* @brief STLMesh provides a simple data structure for manipulating a mesh
 * @note currently only supports triangulated surfaces
 * @note templated on index integer size (int, long int, uint, etc.)
*/
template<typename T>
struct STLMesh {
  std::map<T, std::array<double, 3>> verts;
  std::vector<std::array<T, 3>> ien;
  std::vector<int> featureTags;
};

////////////////////////////////////////////////////////////////////////////////
class MeshMapper {
public:

  /* @brief default constructor */
  MeshMapper() = default;

  /* @brief Get a list of pointers to all the faces on the surface of the geometry
   * @param[out] srfFaces: List of pointers to faces on the surface
   * @param[in] msh: the mesh in question
  */
  static std::vector<apf::MeshEntity*> mapSurfaceFaces(apf::Mesh2* msh);

  /* @brief from a vector of elements, get all the vertices in a map
   * @param[out] vertMap: a map of global index to vertex node
   * @param[in]  list: a list of MeshEntity pointers (i.e. all tris on boundry)
   * @param[in]  gn: a global numbering for the mesh where list comes from
  */
  static std::map<int, std::array<double,3>> mapVertsFromList(
    std::vector<apf::MeshEntity*>& list,
    apf::Numbering* gn
  );

  /* @brief collect list of model features of given topological dimension
   * @param[in] msh = pointer to mesh object
   * @param[in] featureDim = topological dimension of feature being asked for
  */
  static std::vector<int> getModelFeatures(apf::Mesh2* msh, int featureDim);


  /* @brief convert a vectorMesh to an STLMesh, removing duplicate nodes
   * @param[in] vm = VectorMesh object (easy to pass around)
   * @param[out] stlm = STLMesh object (easy to manipulate)
   */
   template <typename T>
   static void mapVectorMeshToSTLMesh(const VectorMesh<T>& vm, STLMesh<T>& stlm);

  /* @brief map the surface of an apf::Mesh2 to an STLMesh
   * @param[in] msh = pointer to mesh object
   * @param[out] stlm = STL mesh
  */
  template <typename T>
  static void mapSurfaceToSTLMesh(apf::Mesh2* msh, STLMesh<T>& stlm);

  /* @brief copy surface mesh data to vector mesh format
   * @param[in] msh = pointer to mesh object
   * @param[out] vm = VectorMesh that is easy to pass around over mpi or network
   */
  template <typename T>
  static void mapSurfaceToVectorMesh(apf::Mesh2* msh, VectorMesh<T>& vm);

  //////////////////////////////////////////////////////////////////////////////
  /* @brief consolidate a VectorMesh onto a single process
   * @note the consolidated mesh will contain duplicates of nodes for every node
   *       that is shared by processes. This is useful, an stlMesh will remove
   *       duplicates, but maintain numbering, and we can use outMsh.vertIndices
   *       along with vertsPerProc as a map to scatter data back to the
   *       processes in the communicator
   * @param[in] inMsh = the parts of the mesh on each process to gather
   * @param[in/out] outMsh = the gathered mesh
                             (only needs to be valid on recvRank process)
   * @param[in/out] vertsPerProc = vector containing nVerts for every process in
                                   communicator
   * @param[in] recvRank = the process rank to receive the mesh
   * @param[in] comm = the MPI communicator to do communication on
  */
  template <typename T>
  static void gatherVectorMesh(
    const VectorMesh<T>& inMsh,
    VectorMesh<T>& outMsh,
    std::vector<int>& vertsPerProc,
    const int recvRank,
    const MPI_Comm comm
  );
private:
};

} // end namespace Shoreline
