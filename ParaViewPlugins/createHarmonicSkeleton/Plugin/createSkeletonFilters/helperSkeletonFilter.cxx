#include "vtkcreateSkeletonFilter.h"
#include <igl/embree/bone_heat.h>

/* There is a nasty compiler error if you include igl/embree/bone_heat.h with
   all the CGAL stuff so this is getting moved to timeout in its own compiled
   file
*/
bool runIGLHeatWeights(
  const Eigen::MatrixXd V,
  const Eigen::MatrixXi F,
  const Eigen::MatrixXd C,
  const Eigen::VectorXi P,
  const Eigen::MatrixXi BE,
  Eigen::MatrixXd& W)
{
  return igl::embree::bone_heat(V, F, C, P, BE, Eigen::MatrixXi(), W);
}
