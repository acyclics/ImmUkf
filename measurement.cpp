#include "measurement.h"
#include "utils.h"

measurement_predict::measurement_predict() {
}

void measurement_predict::initialize(int NSIGMA, double W, double W0_m, double W0_c, ALL_VAR) {
	
	_NSIGMA = NSIGMA;
	_R = MatrixXd::Zero(NZ, NZ);
	_R(0, 0) = VAR_PX;
	_R(1, 1) = VAR_PY;
	_R(2, 2) = VAR_PZ;
	_R(3, 3) = VAR_PROLL;
	_R(4, 4) = VAR_PPITCH;
	_R(5, 5) = VAR_PYAW;

	_WEIGHTS_m.resize(NSIGMA);
	_WEIGHTS_m[0] = W0_m;
	_WEIGHTS_c.resize(NSIGMA);
	_WEIGHTS_c[0] = W0_c;

	for (int i = 1; i < NSIGMA; ++i) {
		_WEIGHTS_m[i] = W;
		_WEIGHTS_c[i] = W;
	}
}

MatrixXd measurement_predict::compute_sigma_z(const MatrixXd& sigma_x) {

	MatrixXd sigma = MatrixXd::Zero(_nz, _NSIGMA);

	for (int c = 0; c < _NSIGMA; ++c) {
		sigma(0, c) = sigma_x(0, c);			// px
		sigma(1, c) = sigma_x(1, c);			// py	
		sigma(2, c) = sigma_x(2, c);			// pz
		sigma(3, c) = sigma_x(3, c);			// proll
		sigma(4, c) = sigma_x(4, c);			// ppitch
		sigma(5, c) = sigma_x(5, c);			// pyaw
	}

	return sigma;
}

MatrixXd measurement_predict::compute_z(const MatrixXd& sigma) {

	VectorXd z = VectorXd::Zero(_nz);

	for (int c = 0; c < _NSIGMA; ++c) {
		z += _WEIGHTS_m[c] * sigma.col(c);
	}

	return z;
}

MatrixXd measurement_predict::compute_S(const MatrixXd& sigma, const MatrixXd& z) {

	VectorXd dz;
	MatrixXd S = MatrixXd::Zero(_nz, _nz);

	for (int c = 0; c < _NSIGMA; ++c) {
		dz = sigma.col(c) - z;
		S += _WEIGHTS_c[c] * dz * dz.transpose();
	}

	S += _R;
	return S;
}

void measurement_predict::process(const MatrixXd& sigma_x) {

	_sigma_z = compute_sigma_z(sigma_x);
	_z = compute_z(_sigma_z);
	_S = compute_S(_sigma_z, _z);
}

VectorXd measurement_predict::getz() const {

	return _z;
}

MatrixXd measurement_predict::getS() const {

	return _S;
}

MatrixXd measurement_predict::get_sigma() const {

	return _sigma_z;
}
