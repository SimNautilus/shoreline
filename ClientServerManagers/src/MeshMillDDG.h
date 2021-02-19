/* @file MeshMillDDG.h
 * @brief Header file for the MeshMillDDG class
 * @details MeshMillDDG provides tools for modifying the geometry of a
 *          distributed finite element mesh. Its features include:
 *          * Operate on a pointer to a distributed mesh.
 *          * Solve a DDG Laplace beltrami problem over the domain given the
 *            prescribed boundary modifications to deform the volume of a mesh
 *
 */
#pragma once
#include <MeshMill.h>
#include <apf.h>
namespace Shoreline {
class MeshMillDDG : public MeshMill {
public:
  MeshMillDDG() = default;
  MeshMillDDG(apf::Mesh2* msh) {
    m_msh = msh;
    setupSolverWithMesh();
  }
  /* @brief form and assemble the harmonic deformation system */
  void formAndAssemble(const std::map<long, std::array<double,3>>& BCs) {
    std::array<std::vector<double>, 6> cotVals;

    elementFormation(cotVals);
    elementAssembly(cotVals, BCs);
  }

  /* @brief Forms the inverse mass matrix to scale Laplace Beltrami operator
   * @details The mass matrix is computed as the barycentric area about each
   *          node in the mesh. As such, the matrix is strictly diagonal and
   *          stored as a petsc Vec object. We really only care about the
   *          inverse so this function computes the inverse of the mass matrix.
   *          This is relatively easy to compute since the mass matrix is
   *          diagonal.
   */
  PetscErrorCode formInvMassMatrix();

  PetscErrorCode solveMeshMotion() override;
  /* @brief destroy the PETSc objects used in the solver */
  PetscErrorCode destroySolverObjects() {
    PetscErrorCode ierr = MeshMill::destroySolverObjects();
    if(m_mIn){ierr = VecDestroy(&m_mIn); CHKERRQ(ierr); m_mIn = NULL;}
    for(Vec& v : m_dxi) {if(v){ierr = VecDestroy(&v); CHKERRQ(ierr); v = NULL;}}
    for(Vec& v : m_rsi) {if(v){ierr = VecDestroy(&v); CHKERRQ(ierr); v = NULL;}}
    return ierr;
  }

  /* @brief apply the computed deformation to the mesh */
  PetscErrorCode applyDeformation() override;

  /* @brief scale the deformation by alpha
   * @note DO THIS BEFORE EVER APPLYING DEFORMATION
  */
  PetscErrorCode scaleDeformation(const double alpha) override;

  /* @brief print diagnostics about the PETSc objects */
  PetscErrorCode printPETScObjectInfo() override {
    PetscErrorCode ierr = 0;
    // Get the rank of the process
    int worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    // get diagnostics about residual vector
    std::array<PetscInt, 3> resLength;
    std::array<std::array<PetscReal,3>, 3> resNorms;
    for(int i = 0; i < 3; ++i) {
      if(m_rsi[i] != PETSC_NULL) {
        ierr = VecGetSize(m_rsi[i], &resLength[i]); CHKERRQ(ierr);
        ierr = VecNorm(m_rsi[i], NORM_1, &resNorms[i][0]); CHKERRQ(ierr);
        ierr = VecNorm(m_rsi[i], NORM_2, &resNorms[i][1]); CHKERRQ(ierr);
        ierr = VecNorm(m_rsi[i], NORM_INFINITY, &resNorms[i][2]); CHKERRQ(ierr);
      }
    }

    // get diagnostics about dy vector
    std::array<PetscInt, 3> dxiLength;
    std::array<std::array<PetscReal,3>, 3> dxiNorms;
    for(int i = 0; i < 3; ++i) {
      if(m_dxi[i] != PETSC_NULL) {
        ierr = VecGetSize(m_dxi[i], &dxiLength[i]); CHKERRQ(ierr);
        ierr = VecNorm(m_dxi[i], NORM_1, &dxiNorms[i][0]); CHKERRQ(ierr);
        ierr = VecNorm(m_dxi[i], NORM_2, &dxiNorms[i][1]); CHKERRQ(ierr);
        ierr = VecNorm(m_dxi[i], NORM_INFINITY, &dxiNorms[i][2]); CHKERRQ(ierr);
      }
    }

    // print from process 1
    if(worldRank == 0) {
      for(int i = 0; i < 3; ++i) {
        std::cout << "Residual Vector " << i << ":\n";
        if(m_rsi[i] == PETSC_NULL) {
          std::cout << "RESIDUAL VECTOR NOT ALLOCATED\n";
        }
        else {
          std::cout << "    length    = " << resLength[i] << "\n";
          std::cout << "    |res|^1   = " << resNorms[i][0] << "\n";
          std::cout << "    |res|^2   = " << resNorms[i][1] << "\n";
          std::cout << "    |res|^inf = " << resNorms[i][2] << "\n";
        }
      }

      for(int i = 0; i < 3; ++i) {
        std::cout << "Solution Vector " << i << ":\n";
        if(m_dxi[i] == PETSC_NULL) {
          std::cout << "SOLUTION VECTOR NOT ALLOCATED\n";
        }
        else {
          std::cout << "    length    = " << dxiLength[i] << "\n";
          std::cout << "    |dxi|^1   = " << dxiNorms[i][0] << "\n";
          std::cout << "    |dxi|^2   = " << dxiNorms[i][1] << "\n";
          std::cout << "    |dxi|^inf = " << dxiNorms[i][2] << "\n";
        }
      }
    }
    return ierr;
  }
private:
  /* @brief Setup the MeshMillDDG solver with info contained in an apf::Mesh2 */
  PetscErrorCode setupSolverWithMesh() override;

  /* @given bcs written into m_dxi, modify m_lhs so that the bcs are cast to
   * the righthand side vector
   */
  PetscErrorCode applyBCsToRHS();

  /* @brief Compute cot(dihedral angles) on each element */
  void elementFormation(std::array<std::vector<double>, 6>& cotVals);

  /* @brief assembles the element matrices into a global matrix */
  PetscErrorCode elementAssembly(
    const std::array<std::vector<double>, 6>& cotVals,
    const std::map<long, std::array<double,3>>& BCs);

  /* @brief inverse of the mass matrix, used to scale Laplace Beltrami matrix */
  Vec m_mIn = NULL;
  /* @brief displacement in the x, y, z directions */
  std::array<Vec, 3> m_dxi = {NULL, NULL, NULL};
  /* @brief rhs vectors in the x, y, z directions */
  std::array<Vec, 3> m_rsi = {NULL, NULL, NULL};
};

} // end namespace Shoreline
