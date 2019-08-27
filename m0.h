#ifndef M0_H
#define M0_H

#include <Eigen/Dense>
#include "ukf.h"
#include <vector>

using namespace Eigen;

ukf ukf0(ALL_DATA, double t);
void model0(MatrixXd& predicted_sigma, const MatrixXd& augmented_sigma, int NSIGMA, double dt);

#endif
