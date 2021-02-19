/* @file MeshMillDDG_ElemAssembly.cpp
 * @details This file contains the MeshMillDDG::elementAssembly routine which
 * assembles element cotan values into global laplace beltrami matrix
 */
// Shoreline includes
#include <MeshMillDDG.h>

// stl includes
#include <vector>

// SCOREC includes
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <apfField.h>

namespace Shoreline {
PetscErrorCode MeshMillDDG::elementAssembly(
  const std::array<std::vector<double>, 6>& cotVals,
  const std::map<long, std::array<double,3>>& BCs)
{
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  PetscErrorCode ierr = 0;

  // zero out the system
  ierr = MatZeroEntries(m_lhs); CHKERRQ(ierr);
  for(int i = 0; i < 3; ++i) {
    ierr = VecZeroEntries(m_rsi[i]); CHKERRQ(ierr);
    ierr = VecZeroEntries(m_dxi[i]);  CHKERRQ(ierr);
  }


  // get element level local numbering to iterate thru cotVals
  apf::Numbering* lnElm = m_msh->findNumbering("LocalElement");
  // get nodal global numbering for matrix assembly
  apf::GlobalNumbering* gn = m_msh->findGlobalNumbering("Global");

  // loop over elements
  apf::MeshIterator* it = m_msh->begin(3);
  apf::MeshEntity* e; // pointer to element
  while((e = m_msh->iterate(it))) {
    if(m_msh->isOwned(e)){
      int elmNum = apf::getNumber(lnElm, e, 0, 0);

      // get the edges of the tet
      apf::Downward edges;
      int nEdges = m_msh->getDownward(e, 1, edges);

      // loop over edges
      for(int edge = 0; edge < nEdges; ++edge) {
        // get the nodes on that edge
        apf::Downward nodes;
        m_msh->getDownward(edges[edge], 0, nodes);

        // get the global numbers of the nodes
        int src = apf::getNumber(gn, nodes[0], 0, 0);
        int dst = apf::getNumber(gn, nodes[1], 0, 0);
        PetscInt ind[2] = {src, dst};
        PetscScalar vals[4] = {-cotVals[edge][elmNum], cotVals[edge][elmNum],
                                cotVals[edge][elmNum], -cotVals[edge][elmNum]};
        // fill values into global matrix
        ierr = MatSetValues(m_lhs, 2,ind, 2,ind, vals,ADD_VALUES);CHKERRQ(ierr);
      }
    }
  }
  m_msh->end(it);
  ierr = MatAssemblyBegin(m_lhs, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_lhs, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  // scale the LaplaceBeltrami operator by the inverse of the mass matrix
  ierr = formInvMassMatrix(); CHKERRQ(ierr);
  ierr = MatDiagonalScale(m_lhs, m_mIn, PETSC_NULL); CHKERRQ(ierr);

  // setup boundary conditions
  std::vector<int> bcRows;
  for(auto bc : BCs) {
    PetscInt nI = bc.first;
    ierr = VecSetValue(m_rsi[0], nI, bc.second[0], INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(m_rsi[1], nI, bc.second[1], INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(m_rsi[2], nI, bc.second[2], INSERT_VALUES);CHKERRQ(ierr);
    bcRows.push_back(nI);
  }
  for(auto& v : m_rsi){
    ierr = VecAssemblyBegin(v); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);   CHKERRQ(ierr);
  }
  ierr = MatZeroRows(m_lhs, bcRows.size(),bcRows.data(), 1.0, NULL, NULL);
  CHKERRQ(ierr);
  return ierr;
}

} // end namespace Shoreline
