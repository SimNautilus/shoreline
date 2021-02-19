/* @file MeshMill_ElemAssembly.cpp
 * @details This file contains the MeshMill::elementAssembly routine which
 * assembles element level stiffness matrices into the global matrix
 */
// Shoreline includes
#include <MeshMill.h>

// stl includes
#include <vector>

// SCOREC includes
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <apfField.h>

namespace Shoreline{
/* @brief Assemble element level stiffness matrices into a global matrix
 */
PetscErrorCode MeshMill::elementAssembly(
  std::vector<double>& kStf,
  std::map<long, std::array<double,3>> BCs)
{
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  PetscErrorCode ierr = 0;

  // zero out the system
  ierr = MatZeroEntries(m_lhs); CHKERRQ(ierr);
  ierr = VecZeroEntries(m_res); CHKERRQ(ierr);
  ierr = VecZeroEntries(m_dy);  CHKERRQ(ierr);


  apf::GlobalNumbering* gn = m_msh->findGlobalNumbering("Global");

  // loop over elements
  apf::MeshIterator* it = m_msh->begin(3);
  int elem = 0;
  apf::MeshEntity* e; // pointer to element
  while((e = m_msh->iterate(it))) {
    // get pointers to the vertices on the face
    apf::Downward ev;
    int nVerts = m_msh->getDownward(e, 0, ev);
    // get the index set of the element
    std::vector<PetscInt> inds(3 * nVerts);
    for(int i = 0; i < nVerts; ++i) {
      int localBlock = 3 * i;
      // int basisI = m_vertIndex2Petsc.at(scorecGlobalNum);
      PetscInt globalI = 3 * apf::getNumber(gn, ev[i], 0, 0);
      inds[localBlock    ] = globalI;
      inds[localBlock + 1] = globalI + 1;
      inds[localBlock + 2] = globalI + 2;
    }
    ierr = MatSetValues(m_lhs, 3 * nVerts, inds.data(), 3 * nVerts,inds.data(),
      &kStf[9*nVerts*nVerts*elem], ADD_VALUES); CHKERRQ(ierr);
    elem++;
  }
  m_msh->end(it);
  ierr = MatAssemblyBegin(m_lhs, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_lhs, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  // fill boundary conditions into rhs vector
  std::vector<int> bcRows;
  int count = 0;
  bcRows.reserve(BCs.size());
  for(auto bc : BCs) {
    PetscInt nI = 3 * bc.first; // block index row start of dof for this node
    PetscInt inds[3] = {nI, nI+1, nI+2};
    ierr = VecSetValues(m_dy, 3, inds, bc.second.data(), INSERT_VALUES); CHKERRQ(ierr);
    bcRows.push_back(nI    );
    bcRows.push_back(nI + 1);
    bcRows.push_back(nI + 2);
    ++count;
  }
  ierr = VecAssemblyBegin(m_dy); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(m_dy); CHKERRQ(ierr);
  ierr = MatZeroRowsColumns(m_lhs, bcRows.size(), bcRows.data(),1.0,m_dy,m_res);
  CHKERRQ(ierr);
  return ierr;
}
}// end namespace Shoreline
