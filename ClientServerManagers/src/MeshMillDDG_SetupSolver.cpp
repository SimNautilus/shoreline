/* @file MeshMillDDG_SetupSolver.cpp
 * @details This file contains routines to setup the MeshMillDDG solver system
 *          used to compute mesh deformations
 */

#include <MeshMillDDG.h>
#include <apfNumbering.h>

#include <numeric>
namespace Shoreline {

PetscErrorCode MeshMillDDG::setupSolverWithMesh() {
  // get an error code for checking the PETSc calls
  PetscErrorCode ierr = 0;

  // get the rank of the process
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  int worldSize;
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

  // get the local numbering of vertices
  apf::Numbering* ownedLocalNumbering = m_msh->findNumbering("Local");
  if(ownedLocalNumbering == NULL) {
    std::cout << "no numbering with name 'Local'\n";
    return -1;
  }
  apf::GlobalNumbering* gn = m_msh->findGlobalNumbering("Global");
  if(gn == NULL) {
    std::cout << "no numbering with name 'Global'\n";
  }

  // get size of locally owned block of matrix
  int nodeCount = apf::countOwned(m_msh, 0);
  std::vector<int> globalCounts(worldSize, 0);
  MPI_Allgather(&nodeCount, 1, MPI_INT, globalCounts.data(), 1, MPI_INT,
    MPI_COMM_WORLD);

  int numRowsGlobal = std::accumulate(globalCounts.begin(),
    globalCounts.end(), 0);

  // number of nonzero blocks in the diagonal square region owned by process
  std::vector<PetscInt> nnzOnDiagSection(nodeCount, 0);
  // number of nonzero blocks in the off diagonal region owned by process
  std::vector<PetscInt> nnzOffDiagSection(nodeCount, 0);


  // loop over vertices on process and count adjacencies
  apf::MeshIterator* it = m_msh->begin(0);
  apf::MeshEntity* v;
  while((v = m_msh->iterate(it))) {
    if(m_msh->isOwned(v)) {
      int localIndex = apf::getNumber(ownedLocalNumbering, v, 0, 0);
      int nNeighbors = m_msh->countUpward(v);
      for(int i = 0; i < nNeighbors; ++i) {
        // look at each edge, figure out if the node adjacent to this node
        // along that edge is owned by the process or not
        // I swear there should be an easier way to do this...
        apf::MeshEntity* edge = m_msh->getUpward(v, i);
        apf::Downward verts;

        m_msh->getDownward(edge, 0, verts);
        // get index of the neighboring vertex along the edge
        int other = (verts[0] == v) ? 1 : 0;
        if(m_msh->isOwned(verts[other])){
          nnzOnDiagSection[localIndex]++;
        } else {
          nnzOffDiagSection[localIndex]++;
        }
      }
    }
  }
  m_msh->end(it);

  // setup the matrix lhs
  ierr = MatCreate(PETSC_COMM_WORLD, &m_lhs); CHKERRQ(ierr);
  ierr = MatSetType(m_lhs, MATMPIAIJ); CHKERRQ(ierr);
  ierr = MatSetSizes(m_lhs, nodeCount, nodeCount,
    numRowsGlobal, numRowsGlobal); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(m_lhs, PETSC_DEFAULT,nnzOnDiagSection.data(),
    PETSC_DEFAULT, nnzOffDiagSection.data());
  CHKERRQ(ierr);
  ierr = MatSetOption(m_lhs, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  CHKERRQ(ierr);

  // Allocate vectors for the solution solve
  ierr = MatCreateVecs(m_lhs, &m_dxi[0], &m_rsi[0]); CHKERRQ(ierr);
  ierr = MatCreateVecs(m_lhs, &m_dxi[1], &m_rsi[1]); CHKERRQ(ierr);
  ierr = MatCreateVecs(m_lhs, &m_dxi[2], &m_rsi[2]); CHKERRQ(ierr);

  // Allocate inverse mass matrix
  ierr = MatCreateVecs(m_lhs, PETSC_NULL, &m_mIn); CHKERRQ(ierr);
  //////////////////////////////////////////////////////////////////////////////
  // setup the vectors used in the solve
  // ierr = MatCreateVecs(m_lhs, &m_dy, &m_res); CHKERRQ(ierr);
  // setup the ksp solver we will start with Jacobi iterations and go from there
  ierr = KSPCreate(PETSC_COMM_WORLD, &m_ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(m_ksp,m_lhs,m_lhs); CHKERRQ(ierr);
  ierr = KSPSetType(m_ksp,KSPGMRES); CHKERRQ(ierr);
  ierr = KSPGetPC(m_ksp,&m_pc); CHKERRQ(ierr);
  ierr = PCSetType(m_pc,PCJACOBI); CHKERRQ(ierr);
  // ierr = PCFactorSetShiftType(m_pc,MAT_SHIFT_POSITIVE_DEFINITE); CHKERRQ(ierr);
  ierr = KSPSetTolerances(m_ksp,
    PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_ksp); CHKERRQ(ierr);

  return ierr;
}

}
