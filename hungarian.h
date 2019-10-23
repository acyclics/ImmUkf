#ifndef HUNGARIAN_H
#define HUNGARIAN_H

#include <vector>
#include <utility>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

vector<int> hungarian_solve(MatrixXd cost);

#endif
