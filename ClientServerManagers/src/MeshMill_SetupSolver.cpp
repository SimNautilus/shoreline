/* @file MeshMill_SetupSolver.cpp
 * @details This file contains routines to setup the MeshMill solver system
 *          used to compute mesh deformations
 */

// Shoreline includes
#include <MeshMill.h>

// SCOREC includes
#include <apfMesh2.h>
#include <apfNumbering.h>
// stl includes
#include <vector>
#include <iostream>
#include <numeric>
#include <algorithm>
// petsc includes
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"

namespace Shoreline {

/* @brief Setup the MeshMill solver with info contained in an apf::Mesh2
 */
PetscErrorCode MeshMill::setupSolverWithMesh() {
  // get an error code for checking the PETSc calls
  PetscErrorCode ierr = 0;

  // dealing with 3D linear elasticity, so we make block sizes of 3
  const PetscInt blockSize = 3;
  // Get the rank of the process
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  int worldSize;
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  // get the rankOffset of the process
  int nodeCount = apf::countOwned(m_msh, 0);
  std::vector<int> globalCounts(worldSize, 0);
  MPI_Allgather(&nodeCount, 1, MPI_INT, globalCounts.data(), 1, MPI_INT,
    MPI_COMM_WORLD);

  //////////////////////////////////////////////////////////////////////////////
  // setup the matrix allocation
  // first get the first matrix row owned by the process
  m_rankOffset = std::accumulate(globalCounts.begin(),
    globalCounts.begin()+worldRank,0);
  // get the number of nodes in the mesh globally.
  int nodeCountGlobal = std::accumulate(globalCounts.begin(),
    globalCounts.end(), 0);

  PetscInt numRowsGlobal = nodeCountGlobal * blockSize;

  // look for the local and global numbering. hopefully theyre setup right
  apf::Numbering* ownedLocalNumbering = m_msh->findNumbering("Local");
  if(ownedLocalNumbering == NULL) {
    std::cout << "no numbering with name 'Local'\n";
    return -1;
  }
  apf::GlobalNumbering* gn = m_msh->findGlobalNumbering("Global");
  if(gn == NULL) {
    std::cout << "no numbering with name 'Global'\n";
  }
  // we are looking at blocks that are 3x3, PETSc requires square blocks

  PetscInt numRowsOnProc = blockSize * nodeCount;

  // number of nonzero blocks in the diagonal square region owned by process
  std::vector<PetscInt> nnzOnDiagSection(nodeCount, 0);
  // number of nonzero blocks in the off diagonal region owned by process
  std::vector<PetscInt> nnzOffDiagSection(nodeCount, 0);

  // there seems to be some weirdness in how SCOREC numbers the indices
  // so let's create a map over vertices to make sure we nail the
  // fill the nonzero count into the two vectors
  // loop over elements on process
  apf::MeshIterator* it = m_msh->begin(3);
  apf::MeshEntity* e;
  PetscInt petscIndex = 0;
  while((e = m_msh->iterate(it))) {
    // get the nodes on the element e
    apf::Downward ev;
    int nVerts = m_msh->getDownward(e, 0, ev);
    // loop over the nodes on the element check if their owned by the process
    for(int i = 0; i < nVerts; ++i) {
      // setup scorec-petsc ordering
      long globalNumber = apf::getNumber(gn, ev[i], 0, 0);
      if(m_vertIndex2Petsc.find(globalNumber) == m_vertIndex2Petsc.end()) {
        m_vertIndex2Petsc.insert({globalNumber, petscIndex + m_rankOffset});
        petscIndex++;
      }
      // check that node i lives on this process
      if(m_msh->isOwned(ev[i]) /*apf::isNumbered(ownedLocalNumbering, ev[i], 0, 0)*/) {
        // get the row corresponding to that node since it lives on the proc
        int localIndex = getNumber(ownedLocalNumbering, ev[i], 0, 0);
        /* now check all dofs that are will be integrated against dof i.
         * if the dof lives on the process, itll get added to the onDiag
         * allocation, otherwise, it'll get added to the offDiag allocation
        */
        for(int j = 0; j < nVerts; ++j) {
          if(m_msh->isOwned(ev[i])/*apf::isNumbered(ownedLocalNumbering, ev[j], 0, 0)*/){
            nnzOnDiagSection[localIndex]++;
          } else {
            nnzOffDiagSection[localIndex]++;
          }
        }
      }
    }
  }
  m_msh->end(it);

  std::for_each(nnzOnDiagSection.begin(), nnzOnDiagSection.end(), [](auto &i){i *= 1.2;});
  std::for_each(nnzOffDiagSection.begin(), nnzOffDiagSection.end(), [](auto &i){i *= 1.2;});

  // setup the matrix lhs
  ierr = MatCreate(PETSC_COMM_WORLD, &m_lhs); CHKERRQ(ierr);
  ierr = MatSetType(m_lhs, MATMPIBAIJ); CHKERRQ(ierr);
  ierr = MatSetSizes(m_lhs, numRowsOnProc, numRowsOnProc,
    numRowsGlobal, numRowsGlobal); CHKERRQ(ierr);
  // ierr = MatMPIBAIJSetPreallocation(m_lhs, blockSize, PETSC_DEFAULT,
  //   nnzOnDiagSection.data(), PETSC_DEFAULT, nnzOffDiagSection.data());
  // CHKERRQ(ierr);


  int d_nz = 1.3* *std::max_element(nnzOnDiagSection.begin(), nnzOnDiagSection.end());
  int o_nz = 3* *std::max_element(nnzOffDiagSection.begin(), nnzOffDiagSection.end());
  ierr = MatMPIBAIJSetPreallocation(m_lhs, blockSize, d_nz,
    PETSC_NULL, o_nz, PETSC_NULL); CHKERRQ(ierr);
  ierr = MatSetOption(m_lhs, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRQ(ierr);
  //////////////////////////////////////////////////////////////////////////////
  // setup the vectors used in the solve
  // ierr = MatCreateVecs(m_lhs, &m_dy, &m_res); CHKERRQ(ierr);
  // setup the ksp solver we will start with Jacobi iterations and go from there
  ierr = KSPCreate(PETSC_COMM_WORLD, &m_ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(m_ksp,m_lhs,m_lhs); CHKERRQ(ierr);
  ierr = KSPSetType(m_ksp, KSPGMRES); CHKERRQ(ierr);
  ierr = KSPGetPC(m_ksp,&m_pc); CHKERRQ(ierr);
  ierr = PCSetType(m_pc,PCBJACOBI); CHKERRQ(ierr);
  // ierr = PCFactorSetShiftType(m_pc,MAT_SHIFT_POSITIVE_DEFINITE); CHKERRQ(ierr);
  ierr = KSPSetTolerances(m_ksp,
    PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_ksp); CHKERRQ(ierr);

  //////////////////////////////////////////////////////////////////////////////
  // setup the ghost ordering of the vector (needed for applying displacement)
  apf::MeshIterator* nodeIt = m_msh->begin(0);
  apf::MeshEntity* v;
  PetscInt ghostIndex = 0;
  int count = 0;
  while((v = m_msh->iterate(nodeIt))) {
    // if this node is a ghosted copy, give it a mapped index
    if(!m_msh->isOwned(v)) {
      ++count;
      long globalIndex = getNumber(gn, v, 0, 0);
      // std::cout << "proc " << worldRank << ": " << globalIndex << "\n";
      m_scorec2GhostNode.insert({globalIndex, ghostIndex});
      ghostIndex++;
    }
  }
  m_msh->end(nodeIt);
  /* now that we have a ghost node map, we need to create a ghost format
   * for the vector m_dy. This will streamline creating an apf::Field
   * object which we will need to apply the displacement computed to the
   * mesh.
  */
  int nGhost = m_scorec2GhostNode.size();
  std::vector<PetscInt> ghosts;
  ghosts.reserve(3 * nGhost);
  std::for_each(m_scorec2GhostNode.begin(), m_scorec2GhostNode.end(),
    [&ghosts](auto node){
      ghosts.push_back(3 * node.first);
      ghosts.push_back(3 * node.first + 1);
      ghosts.push_back(3 * node.first + 2);
    }
  );
  ierr = VecCreateGhost(PETSC_COMM_WORLD, numRowsOnProc, numRowsGlobal,
    ghosts.size(), ghosts.data(), &m_dy); CHKERRQ(ierr);

  ierr = VecDuplicate(m_dy, &m_res); CHKERRQ(ierr);
  return ierr;
}

////////////////////////////////////////////////////////////////////////////////
// @brief just display the matrix sparsity to a window
void MeshMill::testSparsityFromMesh() {

  apf::GlobalNumbering* gn = m_msh->findGlobalNumbering("Global");
  apf::Numbering* ownedLocalNumbering = m_msh->findNumbering("Local");
  // fill the matrix with 1's wherever there is a dof that'd get hit
  // loop over elements on process
  apf::MeshIterator* it = m_msh->begin(3);
  apf::MeshEntity* e;
  while((e = m_msh->iterate(it))) {
    // get the nodes on the element e
    apf::Downward ev;
    int nVerts = m_msh->getDownward(e, 0, ev);
    // loop over the nodes on the element check if their owned by the process
    for(int i = 0; i < nVerts; ++i) {
      // check that node i lives on this process
      if(apf::isNumbered(ownedLocalNumbering, ev[i], 0, 0)) {
        // get the row corresponding to that node since it lives on the proc
        int localI = getNumber(gn, ev[i], 0, 0);
        /* now check all dofs that are will be integrated against dof i.
         * if the dof lives on the process, itll get added to the onDiag
         * allocation, otherwise, it'll get added to the offDiag allocation
        */
        for(int j = 0; j < nVerts; ++j) {
          int localJ = getNumber(gn, ev[j], 0, 0);
          MatSetValue(m_lhs, localI, localJ, 1.0, INSERT_VALUES);
        }
      }
    }
  }

  PetscViewer viewer;
  PetscViewerCreate(MPI_COMM_WORLD, &viewer);
  PetscViewerDrawOpen(MPI_COMM_WORLD, NULL , NULL, 10, 10, PETSC_DECIDE, PETSC_DECIDE, &viewer);
  MatView(m_lhs, viewer);
  PetscViewerDestroy(&viewer);
}
}
