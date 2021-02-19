/* @file MeshMill.h
 * @brief Header file for the MeshMill class
 * @details MeshMill provides tools for modifying the geometry of a distributed
 *          finite element mesh. Its features include:
 *          * Operate on a pointer to a distributed mesh.
 *          * Extract a surface mesh where boundary manipulations can be done
 *          * Solve a linear elasticity problem over the domain given the
 *            prescribed boundary modifications
 *
 */
// #include <apfMesh2.h>
#pragma once
// PETSc includes
#include "petsc.h"
#include "petscmat.h"

// SCOREC includes
#include <apfMesh2.h>

// stl includes
#include <vector>
#include <array>
#include <map>
#include <iostream>

namespace Shoreline {

class MeshMill {
public:
  MeshMill() = default;
  MeshMill(apf::Mesh2* msh) {
    m_msh = msh;
    setupSolverWithMesh();
  }
  ~MeshMill() = default;
  static void print();

  /* @brief form and assemble the linear elasticity system */
  void formAndAssemble(std::map<long, std::array<double,3>> BCs) {
    std::vector<double> kStf;
    elementFormation(kStf);
    elementAssembly(kStf, BCs);
  }

  /* @brief performs the solve of the system produced in formAndAssemble */
  virtual PetscErrorCode solveMeshMotion();

  /* @deprecated */
  void testSparsityFromMesh();

  // @brief destroy the PETSc objects used in the solver
  PetscErrorCode destroySolverObjects() {
    PetscErrorCode ierr = 0;
    // since all these objects are default set to PETSC_NULL, we can
    // check that they are initialized and if so, destroy them.
    if(m_lhs){ierr = MatDestroy(&m_lhs); CHKERRQ(ierr); m_lhs = NULL;}
    if(m_dy) {ierr = VecDestroy(&m_dy);  CHKERRQ(ierr); m_dy  = NULL;}
    if(m_res){ierr = VecDestroy(&m_res); CHKERRQ(ierr); m_res = NULL;}
    if(m_ksp){ierr = KSPDestroy(&m_ksp); CHKERRQ(ierr); m_ksp = NULL;}
    return ierr;
  }

  /* @brief print diagnostics about the PETSc objects */
  virtual PetscErrorCode printPETScObjectInfo() {
    PetscErrorCode ierr = 0;
    // Get the rank of the process
    int worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    // get diagnostics about residual vector
    PetscInt resLength = 0;
    std::array<PetscReal,3> resNorms;

    if(m_res != PETSC_NULL) {
      ierr = VecGetSize(m_res, &resLength); CHKERRQ(ierr);
      ierr = VecNorm(m_res, NORM_1, &resNorms[0]); CHKERRQ(ierr);
      ierr = VecNorm(m_res, NORM_2, &resNorms[1]); CHKERRQ(ierr);
      ierr = VecNorm(m_res, NORM_INFINITY, &resNorms[2]); CHKERRQ(ierr);
    }

    // get diagnostics about dy vector
    PetscInt dyLength = 0;
    std::array<PetscReal,3> dyNorms;
    if(m_dy != PETSC_NULL) {
      ierr = VecGetSize(m_dy, &dyLength); CHKERRQ(ierr);
      ierr = VecNorm(m_dy, NORM_1, &dyNorms[0]); CHKERRQ(ierr);
      ierr = VecNorm(m_dy, NORM_2, &dyNorms[1]); CHKERRQ(ierr);
      ierr = VecNorm(m_dy, NORM_INFINITY, &dyNorms[2]); CHKERRQ(ierr);
    }

    if(worldRank == 0) {
      std::cout << "Residual Vector:\n";
      if(m_res == PETSC_NULL) {
        std::cout << "RESIDUAL VECTOR NOT ALLOCATED\n";
      }
      else {
        std::cout << "    length    = " << resLength << "\n";
        std::cout << "    |res|^1   = " << resNorms[0] << "\n";
        std::cout << "    |res|^2   = " << resNorms[1] << "\n";
        std::cout << "    |res|^inf = " << resNorms[2] << "\n";
      }
      std::cout << "dy Vector:\n";
      if(m_dy == PETSC_NULL) {
        std::cout << "DY VECTOR NOT ALLOCATED\n";
      }
      else {
        std::cout << "    length    = " << dyLength << "\n";
        std::cout << "    |dy|^1   = " << dyNorms[0] << "\n";
        std::cout << "    |dy|^2   = " << dyNorms[1] << "\n";
        std::cout << "    |dy|^inf = " << dyNorms[2] << "\n";
      }
    }
    return ierr;
  }
  /* @bief print the matrix to a matlab .m file */
  PetscErrorCode viewMatrix();

  /* @bief print the solution vector to a matlab .m file */
  PetscErrorCode viewSolution();

  /* @brief apply the computed deformation to the mesh */
  virtual PetscErrorCode applyDeformation();

  /* @brief scale the deformation by alpha
   * @note DO THIS BEFORE EVER APPLYING DEFORMATION
  */
  virtual PetscErrorCode scaleDeformation(const double alpha);
  
protected:
  /* @brief Setup the MeshMill solver with info contained in an apf::Mesh2 */
  virtual PetscErrorCode setupSolverWithMesh();

  /* @brief Compute element level stiffness matrices */
  void elementFormation(std::vector<double>& kStf);

  /* @brief assembles the element stiffness matrices into a matrix */
  PetscErrorCode elementAssembly(std::vector<double>& kStf,
    std::map<long, std::array<double,3>> BCs
  );

  apf::Mesh2* m_msh;
  std::map<long, PetscInt> m_vertIndex2Petsc;

  /* @brief maps a scorec global vertex index to the ghost index of a
   *        localized PETSc vector
   * @details we need a map between global scorec node index to ghost index
   *          because we have no runtime guaruntee that ghosted nodes will be
   *          in order or contiguous, but when one wants the local form of a
   *          PETSc vector, the returned form will be have the ghosted data
   *          contiguous in the array. Here is an excerpt from
   *          VecGhostGetLocalForm:
   *
   *          To update the ghost values from the locations on the other
   *          processes one must call VecGhostUpdateBegin() and
   *          VecGhostUpdateEnd() before accessing the ghost values.
   *          Thus normal usage is
   *
   *          VecGhostUpdateBegin(x,INSERT_VALUES,SCATTER_FORWARD);
   *          VecGhostUpdateEnd(x,INSERT_VALUES,SCATTER_FORWARD);
   *          VecGhostGetLocalForm(x,&xlocal);
   *          VecGetArray(xlocal,&xvalues);
   *           access the non-ghost values in locations xvalues[0:n-1]
   *           and ghost values in locations xvalues[n:n+nghost];
   *          VecRestoreArray(xlocal,&xvalues);
   *          VecGhostRestoreLocalForm(x,&xlocal);
   */
  std::map<long, PetscInt> m_scorec2GhostNode;

  /* @brief the 1st row index in the PETSc matrix region owned by
   * this process(rank)
   */
  int m_rankOffset = 0;

  /* @brief lhs matrix for mesh motion computation */
  Mat m_lhs = PETSC_NULL;

  /* @brief the ksp solver */
  KSP m_ksp = PETSC_NULL;

  /* @brief the preconditioner object */
  PC  m_pc  = PETSC_NULL;

  /* @brief the displacement vector to compute */
  Vec m_dy  = PETSC_NULL;

  /* @brief the solver residual (the rhs vector) */
  Vec m_res = PETSC_NULL;
};

} // end namespace Shoreline
