/* @file MeshMill_ElemFormation.cpp
 * @details This file contains the MeshMill::elementFormation routine
 *          which computes all element stiffness matrices which can be assembled
 *          into a matrix
 */
// Shoreline includes
#include <MeshMill.h>

// stl includes
#include <array>
#include <vector>
#include <iostream>

// SCOREC includes
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <apfField.h>
// mkl includes
#include <mkl_cblas.h>

namespace Shoreline{
/* @brief Compute element level stiffness matrices
 * @note the le_ prefix denotes "linear elasticity"
 *       and refers to material properties
 */
void MeshMill::elementFormation(std::vector<double>& kStf) {
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  // get integration order of accuracy
  const int integrationOrder = 2;
  double kLoc[9];
  // linear elasticity constants
  constexpr std::array<std::array<double,3>, 3> unitVecs = {{
    {{1.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}}, {{0.0, 0.0, 1.0}}
  }};
  // elastic modulus constant
  const double le_eC = 1.0;
  // poisson ratio
  const double le_nu = 0.2;

  // get some properties of the mesh
  int nElems = apf::countOwned(m_msh, 3);
  // std::cout << "nElems on proc = " << nElems << "\n";
  int nVerts = 4;
  // make sure the element stiffness matrix storage is sized properly
  // this is a no-op if it is properly sized already
  kStf.resize(nElems * 9 * nVerts * nVerts);

  // loop over elements
  apf::MeshIterator* it = m_msh->begin(3);
  int elem = 0;
  apf::MeshEntity* e; // pointer to element
  while( (e = m_msh->iterate(it)) ) {
    if(!m_msh->isOwned(e)) continue;

    int elementOffset = 9 * nVerts * nVerts * elem;
    // need a apf::meshElement to get geometry info from the element
    apf::MeshElement* mshElm = apf::createMeshElement(m_msh, e);
    int nGauss = apf::countIntPoints(mshElm, integrationOrder);
    apf::EntityShape* shape = m_msh->getShape()->getEntityShape(m_msh->getType(e));
    
    // std::cout << "entityshape created\n";
    // loop through integration points

    for(int q = 0; q < nGauss; ++q) {
      // get the integration information
      apf::Vector3 quadPt;
      apf::getIntPoint(mshElm, integrationOrder, q, quadPt);
      double quadWt = apf::getIntWeight(mshElm, integrationOrder, q);

      // compute shape function at quadrature point q
      // apf::NewArray<double> elementBases;
      // apf::getValues(m_msh, e, quadPt, elementBases);

      // compute parametric gradients at quad point q
      apf::NewArray<apf::Vector3> localGradients, elmBGrad;
      shape->getLocalGradients(m_msh, e, quadPt, localGradients);
      elmBGrad.allocate(nVerts);

      // compute jacobian at quad point q
      apf::Matrix3x3 jacobian;
      apf::getJacobian(mshElm, quadPt, jacobian);
      apf::Matrix3x3 jinv;
      apf::getJacobianInv(mshElm, quadPt, jinv);
      const double detJ  = apf::getJacobianDeterminant(jacobian, 3);
      const double wDetJ = quadWt * detJ;
      
      // compute gradients in physical space
      for (int i=0; i < nVerts; ++i)
        elmBGrad[i] = jinv * localGradients[i];

      // element elastic modulus is proportional to 1/sqrt of jacobian
      const double le_E = le_eC * (1.0 / sqrt(detJ));

      // element material properties
      const double le_Lambda = le_E / ((1.0 + le_nu) * (1.0 - 2.0 * le_nu));
      const double le_LxNu = le_Lambda * le_nu;
      const double le_Lx1m2nu = le_Lambda * (1.0 - 2.0 * le_nu);

      // compute linear elastic modulus matrix, D
      const double le_D[36] = {
        le_Lambda * (1 - le_nu), le_LxNu,                 le_LxNu,                 0.0,        0.0,        0.0,
        le_LxNu,                 le_Lambda * (1 - le_nu), le_LxNu,                 0.0,        0.0,        0.0,
        le_LxNu,                 le_LxNu,                 le_Lambda * (1 - le_nu), 0.0,        0.0,        0.0,
        0.0,                     0.0,                     0.0,                     le_Lx1m2nu, 0.0,        0.0,
        0.0,                     0.0,                     0.0,                     0.0,        le_Lx1m2nu, 0.0,
        0.0,                     0.0,                     0.0,                     0.0,        0.0,        le_Lx1m2nu
      };

      // loop over shape functions in i
      for(int i = 0; i < nVerts; ++i) {
        const double b_aTranspose[18] = {
          elmBGrad[i][0], 0.0, 0.0, elmBGrad[i][1], 0.0, elmBGrad[i][2],
          0.0, elmBGrad[i][1], 0.0, elmBGrad[i][0], elmBGrad[i][2], 0.0,
          0.0, 0.0, elmBGrad[i][2], 0.0, elmBGrad[i][1], elmBGrad[i][0]
        };

        // compute b_a^T * D * b_b * qwt * detJ
        // 1) compute tmp = b_a^T * D = 3 x 6
        double tmp[18];
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
          3, 6, 6, 1.0, b_aTranspose, 6, le_D, 6, 0.0, tmp, 6);

        // loop over shape functions in j
        for(int j = 0; j < nVerts; ++j) {
          const double b_b[18] = {
            elmBGrad[j][0], 0.0,            0.0,
            0.0,            elmBGrad[j][1], 0.0,
            0.0,            0.0,            elmBGrad[j][2],
            elmBGrad[j][1], elmBGrad[j][0], 0.0,
            0.0,            elmBGrad[j][2], elmBGrad[j][1],
            elmBGrad[j][2], 0.0,            elmBGrad[j][0]
          };

          // 3) compute kLoc = kLoc + tmp*b_b * wDetJ
          cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            3, 3, 6, wDetJ, tmp, 6, b_b, 3, 0.0, kLoc, 3);

          // loop over directions and combine into element stiffness matrix
          for(int di = 0; di < 3; ++di) {
            int a = 3 * i + di;
            for(int dj = 0; dj < 3; ++dj) {
              int b = 3 * j + dj;
              int eI = (elementOffset) + (a + b * 3 * nVerts);
              kStf[eI] += unitVecs[dj][0] * (kLoc[0] * unitVecs[di][0] +
                                             kLoc[1] * unitVecs[di][1] +
                                             kLoc[2] * unitVecs[di][2] )
                        + unitVecs[dj][1] * (kLoc[3] * unitVecs[di][0] +
                                             kLoc[4] * unitVecs[di][1] +
                                             kLoc[5] * unitVecs[di][2] )
                        + unitVecs[dj][2] * (kLoc[6] * unitVecs[di][0] +
                                             kLoc[7] * unitVecs[di][1] +
                                             kLoc[8] * unitVecs[di][2] );
            }
          }
        }
      }
    }
    elem++;
    apf::destroyMeshElement(mshElm);
  }
}
}// end namespace Shoreline
