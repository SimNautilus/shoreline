#include "MeshMapper.h"
#include "ShorelineMPI.h"
#include <iostream>
#include <apfNumberingClass.h>

namespace Shoreline{
// #include <apfVectorField.h>
////////////////////////////////////////////////////////////////////////////////
std::vector<apf::MeshEntity*> MeshMapper::mapSurfaceFaces(apf::Mesh2* msh) {
  std::vector<apf::MeshEntity*> srfFaces;
  // get list of model surfaces
  std::vector<int> modelSrfs = {};
  gmi_model* model  = msh->getModel();
  gmi_iter*  geoItr;
  gmi_ent*   entity;
  geoItr = gmi_begin(model, 2); // 2 for surfaces
  while( (entity = gmi_next(model, geoItr)) ) {
    modelSrfs.push_back(gmi_tag(model, entity));
  }
  gmi_end(model, geoItr); // free the iterator

  // loop over faces on partition and get a list of those on the model surfaces
  apf::MeshIterator* it = msh->begin(2);
  apf::MeshEntity* f;
  while( (f = msh->iterate(it)) ) {
    int fTag = msh->getModelTag(msh->toModel(f));
    // check if the element tag lands on the geometric surface
    if(modelSrfs.end() != std::find(modelSrfs.begin(), modelSrfs.end(), fTag)) {
      srfFaces.push_back(f);
    }
  }
  msh->end(it); // free the iterator

  // return the list of faces
  return srfFaces;
}

////////////////////////////////////////////////////////////////////////////////
std::map<int, std::array<double,3>> MeshMapper::mapVertsFromList(
  std::vector<apf::MeshEntity*>& list, apf::Numbering* gn)
{
  apf::Mesh* msh = gn->getMesh();
  std::map<int,std::array<double,3>> vertMap; // vertex map to return

  // get a global nodes list
  apf::DynamicArray<apf::Node> nodes;
  apf::getNodes(gn, nodes);

  // get the nodal indices of the entities in the list
  for(auto& f : list) {
    int nen = apf::getShape(gn)->getEntityShape(msh->getType(f))->countNodes();
    std::cout << "nen = " << nen << "\n\n";

    // get the vertices and put them in the map
    std::array<double,3> vertData; // assume triangle faces
    apf::Vector3 vertRaw;
    for(int i = 0; i < nen; ++i) {
      int ind = gn->get(f, i, 0);
      std::cout << ind << "\n";
      msh->getPoint(nodes[ind].entity,0,vertRaw);
      std::copy(&vertRaw[0],&vertRaw[0]+3, &vertData[0]);
      vertMap.insert({ind, vertData});
    }
  }

  return vertMap;
}
////////////////////////////////////////////////////////////////////////////////
std::vector<int> MeshMapper::getModelFeatures(apf::Mesh2* msh, int featureDim) {
  std::vector<int> modelFeatures;
  gmi_model* model  = msh->getModel();
  gmi_iter*  geoItr;
  gmi_ent*   entity;
  geoItr = gmi_begin(model, featureDim); // 2 for surfaces
  while( (entity = gmi_next(model, geoItr)) ) {
    modelFeatures.push_back(gmi_tag(model, entity));
  }
  gmi_end(model, geoItr); // free the iterator
  return modelFeatures;
}

////////////////////////////////////////////////////////////////////////////////
template<typename T>
void MeshMapper::mapSurfaceToSTLMesh(apf::Mesh2* msh, STLMesh<T>& stlm) {
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  stlm.ien.resize(0);

  apf::GlobalNumbering* gn = msh->findGlobalNumbering("Global");
  // get list of model surfaces
  std::vector<int> modelSrfs = getModelFeatures(msh , 2);

  // loop over faces on partition and get a list of those on the model surfaces
  apf::MeshIterator* it = msh->begin(2);
  apf::MeshEntity* f;

  // get the vertices and put them in the map
  std::array<double,3> vertData; // assume triangle faces
  apf::Vector3 vertRaw;
  std::array<T,3> faceInds;
  int ownedSrfNodes = 0;
  // loop over all faces in partition
  while( (f = msh->iterate(it)) ) {
    int fTag = msh->getModelTag(msh->toModel(f));
    // check if the element tag lands on the geometric surface
    if(modelSrfs.end() != std::find(modelSrfs.begin(), modelSrfs.end(), fTag)) {
      // get pointers to the vertices on the face
      apf::Downward fv;
      int nVerts = msh->getDownward(f, 0, fv);
      // loop over vertices
      for(int i = 0; i < nVerts; ++i) {
        // get global index of the vertex
        long ind = getNumber(gn, fv[i], 0, 0);
        msh->getPoint(fv[i], 0, vertRaw);
        if(!stlm.verts.count(ind) && msh->isOwned(fv[i])) {
          ownedSrfNodes++;
        }
        // get the data for the vertex, copy into stl array
        std::copy(&vertRaw[0],&vertRaw[0]+3, &vertData[0]);
        // insert the vertex data into the std::map
        stlm.verts.insert({ind, vertData});

        faceInds[i] = ind;
      }
      stlm.ien.push_back(faceInds);
      stlm.featureTags.push_back(fTag);
    }
  }
  msh->end(it); // free the iterator
}


template
void MeshMapper::mapSurfaceToSTLMesh<int>(apf::Mesh2* msh, STLMesh<int>& stlm);
template
void MeshMapper::mapSurfaceToSTLMesh<long>(apf::Mesh2* msh, STLMesh<long>& stlm);
////////////////////////////////////////////////////////////////////////////////
template <typename T>
void MeshMapper::mapSurfaceToVectorMesh(apf::Mesh2* msh, VectorMesh<T>& vm) {
  // first make an stl mesh so we don't record a bunch of duplicate values
  STLMesh<T> stlm;
  mapSurfaceToSTLMesh(msh, stlm);

  // record the number of vertices and faces on the surface mesh
  vm.nVerts = stlm.verts.size();
  vm.nFaces = stlm.ien.size();

  // extract vertex indices from map into a contiguous array
  // extract vertex coordinates from map into a contiguous array
  vm.vertIndices.resize(0);
  vm.vertCoords.resize(0);
  for(auto p : stlm.verts) {
    vm.vertIndices.push_back(p.first);
    vm.vertCoords.push_back(p.second[0]);
    vm.vertCoords.push_back(p.second[1]);
    vm.vertCoords.push_back(p.second[2]);
  }

  // extract feature tags from stl mesh
  vm.featureTags = stlm.featureTags;

  // extract face indices from ien into contiguous array
  vm.faceIndices.resize(0);
  for(auto p : stlm.ien) {
    vm.faceIndices.push_back(p[0]);
    vm.faceIndices.push_back(p[1]);
    vm.faceIndices.push_back(p[2]);
  }
}

// template specializations
template void MeshMapper::mapSurfaceToVectorMesh<long>(
  apf::Mesh2* msh, VectorMesh<long>& vm
);
template void MeshMapper::mapSurfaceToVectorMesh<int>(
  apf::Mesh2* msh, VectorMesh<int>& vm
);

////////////////////////////////////////////////////////////////////////////////
template <typename T>
void MeshMapper::mapVectorMeshToSTLMesh(
  const VectorMesh<T>& vm, STLMesh<T>& stlm)
{
  for(int i = 0; i < vm.nVerts; ++i) {
    std::array<double, 3> coord =
      {vm.vertCoords[3*i], vm.vertCoords[3*i+1], vm.vertCoords[3*i+2]};
    stlm.verts.insert({vm.vertIndices[i], coord});
  }

  stlm.featureTags.resize(0);
  stlm.featureTags.reserve(vm.nFaces);
  for(int i = 0; i < vm.nFaces; ++i) {
    std::array<T, 3> face =
      {vm.faceIndices[3*i], vm.faceIndices[3*i+1], vm.faceIndices[3*i+2]};
    stlm.ien.push_back(face);
    stlm.featureTags.push_back(vm.featureTags[i]);
  }
  // stlm.featureTags = vm.featureTags;
}

// template specializations
template void MeshMapper::mapVectorMeshToSTLMesh<long>(
  const VectorMesh<long>& vm, STLMesh<long>& stlm);
template void MeshMapper::mapVectorMeshToSTLMesh<int>(
  const VectorMesh<int>& vm, STLMesh<int>& stlm);

////////////////////////////////////////////////////////////////////////////////
template <typename T>
void MeshMapper::gatherVectorMesh(
  const VectorMesh<T>& inMsh,
  VectorMesh<T>& outMsh,
  std::vector<int>& vertsPerProc,
  const int recvRank,
  const MPI_Comm comm)
{
  // get how big the communicator network is
  int commSize;
  MPI_Comm_size(comm, &commSize);
  // Get the rank of the process
  int commRank;
  MPI_Comm_rank(comm, &commRank);

  // get the index datatype so MPI can know how much data to pass
  MPI_Datatype mpi_type = MPI_LONG;
  if(std::is_same<int, T>::value) { mpi_type = MPI_INT; }

  // get how many elements are getting sent per process
  // (useful for mapping vectors back to processes)
  vertsPerProc.resize((commRank==recvRank) ? commSize : 0);
  MPI_Gather(&inMsh.nVerts, 1, MPI_INT,
    vertsPerProc.data(), 1, MPI_INT, recvRank, comm);

  // gather all the vertices, vertex indices and face connectivity on proc 0
  Shoreline::gatherVectors(inMsh.vertIndices,outMsh.vertIndices,recvRank,comm);
  Shoreline::gatherVectors(inMsh.vertCoords, outMsh.vertCoords, recvRank,comm);
  Shoreline::gatherVectors(inMsh.faceIndices,outMsh.faceIndices,recvRank,comm);
  Shoreline::gatherVectors(inMsh.featureTags,outMsh.featureTags,recvRank,comm);

  outMsh.nVerts = outMsh.vertIndices.size();
  outMsh.nFaces = outMsh.faceIndices.size()/3;
}
// template specializations of the gatherVectorMesh function
template void MeshMapper::gatherVectorMesh<int>(
  const VectorMesh<int>& inMsh,
  VectorMesh<int>& outMsh,
  std::vector<int>& vertsPerProc,
  const int recvRank,
  const MPI_Comm comm);
template void MeshMapper::gatherVectorMesh<long>(
  const VectorMesh<long>& inMsh,
  VectorMesh<long>& outMsh,
  std::vector<int>& vertsPerProc,
  const int recvRank,
  const MPI_Comm comm);
// void MeshMapper::subMeshFromList(const std::vector<apf::MeshEntity*> list, )
} // end namespace Shoreline
