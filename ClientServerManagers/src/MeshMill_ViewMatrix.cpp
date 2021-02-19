/* @file MeshMill_ViewMatrix.cpp
 * @brief simple viewer code to inspect the matrix m_lhs in MeshMill
 */

#include <MeshMill.h>

#include "petscdraw.h"
#include "petscviewer.h"

#include <iostream>

namespace Shoreline {
PetscErrorCode MeshMill::viewMatrix() {
  PetscErrorCode ierr = 0;
  PetscViewer viewer;
  ierr = PetscObjectSetName((PetscObject)m_lhs, "m_lhs"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "../output/lhs.m", &viewer);CHKERRQ(ierr);

  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = MatView(m_lhs,viewer);CHKERRQ(ierr);
  return ierr;
}

PetscErrorCode MeshMill::viewSolution() {
  PetscErrorCode ierr = 0;
  PetscViewer viewer;
  ierr = PetscObjectSetName((PetscObject)m_dy, "m_dy"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "../output/dy.m", &viewer);CHKERRQ(ierr);

  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = VecView(m_dy,viewer);CHKERRQ(ierr);
  return ierr;
}
}
