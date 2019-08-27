#ifndef M1_H
#define M1_H

#include <Eigen/Dense>
#include "ukf.h"
#include <vector>

using namespace Eigen;

ukf ukf1(ALL_DATA, double t);
void model1(MatrixXd& predicted_sigma, const MatrixXd& augmented_sigma, int NSIGMA, double dt);

#endif
