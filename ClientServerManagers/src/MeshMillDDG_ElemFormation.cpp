/* @file MeshMillDDG_ElemFormation.cpp
 * @brief Contains functions related to local element matrix formation
 *        for the MeshMillDDG mesh deformation solver
 */
// Shoreline includes
#include <MeshMillDDG.h>

// mkl includes
#include <mkl.h>

// stl includes
#include <array>
#include <vector>
#include <algorithm>

// SCOREC includes
#include <apfNumbering.h>
#include <apf.h>

namespace Shoreline {

/* @brief return the square of the input */
template <typename T> inline T sqr(T a) {return a * a;}

void computeDihedralAngles(
  const std::array<std::vector<double>, 6> eL, // edge lengths
  const std::array<std::vector<double>, 4> fA, // face areas
  std::array<std::vector<double>, 6>& cosTheta // cos(dihedral angles)
)
{
  const int numel = eL[0].size();
  double Hsqr[6];

  // loop over elements
  for(int e = 0; e < numel; ++e) {
    // compute Hsqr
    Hsqr[0] = 0.0625 * (4. * sqr(eL[3][e]) * sqr(eL[0][e])
                       - sqr(sqr(eL[1][e]) + sqr(eL[4][e])
                           - sqr(eL[2][e]) - sqr(eL[5][e])));
    Hsqr[1] = 0.0625 * (4. * sqr(eL[4][e]) * sqr(eL[1][e])
                       - sqr(sqr(eL[2][e]) + sqr(eL[5][e])
                           - sqr(eL[3][e]) - sqr(eL[0][e])));
    Hsqr[2] = 0.0625 * (4. * sqr(eL[5][e]) * sqr(eL[2][e])
                       - sqr(sqr(eL[3][e]) + sqr(eL[0][e])
                           - sqr(eL[4][e]) - sqr(eL[1][e])));
    Hsqr[3] = 0.0625 * (4. * sqr(eL[0][e]) * sqr(eL[3][e])
                       - sqr(sqr(eL[4][e]) + sqr(eL[1][e])
                           - sqr(eL[5][e]) - sqr(eL[2][e])));
    Hsqr[4] = 0.0625 * (4. * sqr(eL[1][e]) * sqr(eL[4][e])
                       - sqr(sqr(eL[5][e]) + sqr(eL[2][e])
                           - sqr(eL[0][e]) - sqr(eL[3][e])));
    Hsqr[5] = 0.0625 * (4. * sqr(eL[2][e]) * sqr(eL[5][e])
                       - sqr(sqr(eL[0][e]) + sqr(eL[3][e])
                           - sqr(eL[1][e]) - sqr(eL[4][e])));
    // compute cosTheta
    cosTheta[0][e] = (Hsqr[0] - sqr(fA[1][e]) - sqr(fA[2][e]))
                   / (-2. * fA[1][e] * fA[2][e]);
    cosTheta[1][e] = (Hsqr[1] - sqr(fA[2][e]) - sqr(fA[0][e]))
                   / (-2. * fA[2][e] * fA[0][e]);
    cosTheta[2][e] = (Hsqr[2] - sqr(fA[0][e]) - sqr(fA[1][e]))
                   / (-2. * fA[0][e] * fA[1][e]);
    cosTheta[3][e] = (Hsqr[3] - sqr(fA[3][e]) - sqr(fA[0][e]))
                   / (-2. * fA[3][e] * fA[0][e]);
    cosTheta[4][e] = (Hsqr[4] - sqr(fA[3][e]) - sqr(fA[1][e]))
                   / (-2. * fA[3][e] * fA[1][e]);
    cosTheta[5][e] = (Hsqr[5] - sqr(fA[3][e]) - sqr(fA[2][e]))
                   / (-2. * fA[3][e] * fA[2][e]);
  }
}

////////////////////////////////////////////////////////////////////////////////
void MeshMillDDG::elementFormation(std::array<std::vector<double>, 6>& cotVals)
{
  // get local numbering
  apf::Numbering* ln = m_msh->findNumbering("LocalElement");
  // get number of elements owned by process
  const int numel = apf::countOwned(m_msh, 3);
  // collect edge lengths on each element. tets have 6 edges
  std::array<std::vector<double>, 6> eL;
  std::for_each(eL.begin(), eL.end(), [numel](auto& v){v.resize(numel);});

  // collect face areas on each element. tets have 4 faces
  std::array<std::vector<double>, 4> fA;
  std::for_each(fA.begin(), fA.end(), [numel](auto& v){v.resize(numel);});
  // collect volume of each element
  std::vector<double> volumes(numel);

  // loop over elements
  apf::MeshIterator* it = m_msh->begin(3);
  apf::MeshEntity* elm;
  while((elm = m_msh->iterate(it))) {
    if(m_msh->isOwned(elm)){
      const int e = apf::getNumber(ln, elm, 0, 0);
      volumes[e] = apf::measure(m_msh, elm);
      apf::Downward faces;
      apf::Downward edges;
      int nFaces = m_msh->getDownward(elm, 2, faces);
      int nEdges = m_msh->getDownward(elm, 1, edges);

      for(int face = 0; face < nFaces; ++face) {
        fA[face][e] = apf::measure(m_msh, faces[face]);
      }
      for(int edge = 0; edge < nEdges; ++edge) {
        eL[edge][e] = apf::measure(m_msh, edges[edge]);
      }
    }
  }
  m_msh->end(it);
  // compute dihedral angles for each node on each element
  std::array<std::vector<double>, 6> cosTheta;
  std::for_each(cosTheta.begin(), cosTheta.end(), [numel](auto& v){v.resize(numel);});
  computeDihedralAngles(eL, fA, cosTheta);
  // compute sines of the dihedral angles
  std::array<std::vector<double>, 6> sinTheta;
  std::for_each(sinTheta.begin(), sinTheta.end(), [numel](auto& v){v.resize(numel);});

  for(int e = 0; e < numel; ++e) {
    sinTheta[0][e] = 3. * volumes[e] * eL[0][e] / (2. * fA[1][e] * fA[2][e]);
    sinTheta[1][e] = 3. * volumes[e] * eL[1][e] / (2. * fA[2][e] * fA[0][e]);
    sinTheta[2][e] = 3. * volumes[e] * eL[2][e] / (2. * fA[0][e] * fA[1][e]);
    sinTheta[3][e] = 3. * volumes[e] * eL[3][e] / (2. * fA[3][e] * fA[0][e]);
    sinTheta[4][e] = 3. * volumes[e] * eL[4][e] / (2. * fA[3][e] * fA[1][e]);
    sinTheta[5][e] = 3. * volumes[e] * eL[5][e] / (2. * fA[3][e] * fA[2][e]);
  }
  // compute entries of cotangent matrix per edge scaled by 1/6 edge length
  std::vector<double> scratch(numel);
  for(int i = 0; i < 6; ++i) {
    cotVals[i].resize(numel);
    vdDiv(numel, cosTheta[i].data(), sinTheta[i].data(), scratch.data());
    vdMul(numel, scratch.data(), eL[i].data(), cotVals[i].data());
    std::for_each(cotVals[i].begin(), cotVals[i].end(), [](auto & p){p /= 6.;});
  }
}

////////////////////////////////////////////////////////////////////////////////
PetscErrorCode MeshMillDDG::formInvMassMatrix() {
  PetscErrorCode ierr = 0;

  // zero out the entries in the matrix
  ierr = VecZeroEntries(m_mIn); CHKERRQ(ierr);

  // get the global numbering of the nodes
  apf::GlobalNumbering* gn = m_msh->findGlobalNumbering("Global");

  // loop over nodes
  apf::MeshIterator* it = m_msh->begin(3);
  apf::MeshEntity* elm;
  // we want to compute the barycentric volume about each node
  // this is actually super easy to do. just get the volume of each element
  // and add 1/4 of that to each vertex volume
  while((elm = m_msh->iterate(it))) {
    if(m_msh->isOwned(elm)) {
      apf::Downward verts;
      int nVerts = m_msh->getDownward(elm, 0, verts);
      // get the volume of the element
      double elemMass = apf::measure(m_msh, elm);
      PetscInt inds[4];
      PetscScalar massInvs[4];
      for(int i = 0; i < nVerts; ++i) {
        inds[i] = apf::getNumber(gn, verts[i], 0, 0);
        // we want 1/mass of element
        massInvs[i] = 4.0 / elemMass;
      }
      // insert into inverse mass matrix
      ierr = VecSetValues(m_mIn, 4, inds, massInvs, ADD_VALUES); CHKERRQ(ierr);
    }
  }
  m_msh->end(it);
  ierr = VecAssemblyBegin(m_mIn); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(m_mIn); CHKERRQ(ierr);
  return ierr;
}
} // end namespace Shoreline
