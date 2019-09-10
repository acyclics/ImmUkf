#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "global_settings.h"
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

class measurement_predict {
	int _nz = NZ;
	int _NSIGMA;
	MatrixXd _R;
	VectorXd _z;
	MatrixXd _S;
	MatrixXd _sigma_z;

	vector<double> _WEIGHTS_m, _WEIGHTS_c;

	MatrixXd compute_sigma_z(const MatrixXd& sigma_x);
	MatrixXd compute_z(const MatrixXd& sigma_z);
	MatrixXd compute_S(const MatrixXd& sigma_z, const MatrixXd& predicted_z);

	public:
	measurement_predict();
	void initialize(int NSIGMA, double W, double W0_m, double W0_c, ALL_VAR);
	void process(const MatrixXd& sigma_x);
	MatrixXd get_sigma() const;
	MatrixXd getS() const;
	VectorXd getz() const;
};

#endif
