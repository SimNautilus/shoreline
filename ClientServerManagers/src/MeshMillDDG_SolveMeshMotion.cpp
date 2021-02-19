/* @file MeshMillDDG_SolveMeshMotion.cpp
 * @details This file contains the MeshMillDDG routine for solving for the
 *          deformation of a volume mesh given a mesh and boundary conditions.
 */

// Shoreline includes
#include <MeshMillDDG.h>

//stl includes
#include <iostream>

// SCOREC includes
#include <apfMesh2.h>
#include <apfField.h>
#include <apfFieldData.h>
#include <apfNumbering.h>

namespace Shoreline {

/* @brief Compute the deformation of a 3D mesh given boundary conditions
 * @details This function computes the deformation of a 3D tetrahedral mesh
 *          in parallel given a set of boundary conditions. The passed in mesh
 *          is updated with the deformation so as to not make a copy of a large
 *          mesh.
 */
PetscErrorCode MeshMillDDG::solveMeshMotion() {
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  PetscErrorCode ierr = 0;
  KSPConvergedReason reason;
  PetscInt its = 0;
  for(int i = 0; i < 3; ++i) {
    ierr = KSPSolve(m_ksp, m_rsi[i], m_dxi[i]); CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(m_ksp, &its); CHKERRQ(ierr);

    if(worldRank == 0) {

    }

    ierr = KSPGetConvergedReason(m_ksp, &reason); CHKERRQ(ierr);
    if(worldRank == 0){
      std::cout << "Direction " << i <<  ", Petsc Iterations for Mesh Motion: "
                << its << ". Converged with reason: ";
      switch (reason) {
        case 0:   std::cout << "KSP_CONVERGED_ITERATING\n";       break;
        case 1:   std::cout << "KSP_CONVERGED_RTOL_NORMAL\n";     break;
        case 2:   std::cout << "KSP_CONVERGED_RTOL\n";            break;
        case 3:   std::cout << "KSP_CONVERGED_ATOL\n";            break;
        case 4:   std::cout << "KSP_CONVERGED_ITS\n";             break;
        case 5:   std::cout << "KSP_CONVERGED_CG_NEG_CURVE\n";    break;
        case 6:   std::cout << "KSP_CONVERGED_CG_CONSTRAINED\n";  break;
        case 7:   std::cout << "KSP_CONVERGED_STEP_LENGTH\n";     break;
        case 8:   std::cout << "KSP_CONVERGED_HAPPY_BREAKDOWN\n"; break;
        case 9:   std::cout << "KSP_CONVERGED_ATOL_NORMAL\n";     break;
        case -2:  std::cout << "KSP_DIVERGED_NULL\n";             break;
        case -3:  std::cout << "KSP_DIVERGED_ITS\n";              break;
        case -4:  std::cout << "KSP_DIVERGED_DTOL\n";             break;
        case -5:  std::cout << "KSP_DIVERGED_BREAKDOWN\n";        break;
        case -6:  std::cout << "KSP_DIVERGED_BREAKDOWN_BICG\n";   break;
        case -7:  std::cout << "KSP_DIVERGED_NONSYMMETRIC\n";     break;
        case -8:  std::cout << "KSP_DIVERGED_INDEFINATE_PC\n";    break;
        case -9:  std::cout << "KSP_DIVERGED_NANORINF\n";         break;
        case -10: std::cout << "KSP_DIVERGED_INDEFINITE_MAT\n";   break;
        case -11: std::cout << "KSP_DIVERGED_PCSETUP_FAILED\n";   break;
      }
    }
  }

  return ierr;
}

////////////////////////////////////////////////////////////////////////////////
/* @brief scale the deformation by alpha
 * @note DO THIS BEFORE EVER APPLYING DEFORMATION
*/
PetscErrorCode MeshMillDDG::scaleDeformation(const double alpha) {
  PetscErrorCode ierr = VecScale(m_dxi[0], alpha); CHKERRQ(ierr);
  ierr = VecScale(m_dxi[1], alpha); CHKERRQ(ierr);
  ierr = VecScale(m_dxi[2], alpha); CHKERRQ(ierr);
  return ierr;
}

////////////////////////////////////////////////////////////////////////////////
/* @brief apply the computed deformation from m_dy onto m_msh */
PetscErrorCode MeshMillDDG::applyDeformation() {
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  PetscErrorCode ierr = 0;

  // get an array of local values of the displacement field
  std::array<Vec, 3> local_dxi = {NULL, NULL, NULL};
  std::array<VecScatter, 3> scatters;
  std::array<PetscScalar*, 3> dispArrays;
  for(int i = 0; i < 3; ++i){
    ierr = VecScatterCreateToAll(m_dxi[i], &scatters[i], &local_dxi[i]);
    CHKERRQ(ierr);
    ierr = VecScatterBegin(scatters[i],m_dxi[i], local_dxi[i],
      INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(scatters[i], m_dxi[i], local_dxi[i],
      INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

    // get the data into a useable format
    ierr = VecGetArray(local_dxi[i], &dispArrays[i]); CHKERRQ(ierr);
  }


  // setup an apf::Field to apply displacements to
  apf::Field* d = m_msh->findField("displacement");
  if(!d) d = apf::createFieldOn(m_msh, "displacement", apf::VECTOR);

  // get the global numbering of the mesh
  apf::GlobalNumbering* gn = m_msh->findGlobalNumbering("Global");
  // loop over vertices and apply displacement to the Field object
  apf::MeshIterator* it = m_msh->begin(0);
  apf::MeshEntity* v;
  while( (v = m_msh->iterate(it)) ) {
    long g = apf::getNumber(gn, v, 0, 0);
    double dx = dispArrays[0][g];
    double dy = dispArrays[1][g];
    double dz = dispArrays[2][g];
    apf::Vector3 disp(dx, dy, dz);
    setVector(d, v, 0, disp);
  }
  m_msh->end(it); // free the iterator
  apf::displaceMesh(m_msh, d);
  apf::synchronizeFieldData<double>(m_msh->getCoordinateField()->getData(),NULL);
  m_msh->acceptChanges();
  // restore the vector's guts to the vector
  MPI_Barrier(MPI_COMM_WORLD);
  // apf::destroyField(d);

  // clean up the petsc objects
  for(int i = 0; i < 3; ++i) {
    ierr = VecRestoreArray(local_dxi[i], &dispArrays[i]); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatters[i]); CHKERRQ(ierr);
    ierr = VecDestroy(&local_dxi[i]); CHKERRQ(ierr);
  }

  return ierr;
}
}
