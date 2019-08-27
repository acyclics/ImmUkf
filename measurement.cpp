#include "measurement.h"
#include "utils.h"

measurement_predict::measurement_predict() {
}

void measurement_predict::initialize(int NSIGMA, double W, double W0_m, double W0_c, double VAR_PX, double VAR_PY, double VAR_PZ, double VAR_PYAW, double VAR_PPITCH, double VAR_PROLL) {
	
	_NSIGMA = NSIGMA;
	_R = MatrixXd(NZ, NZ);

	_R << VAR_PX,     0,     0,     0,     0,     0,
		  0 , VAR_PY,     0,     0,     0,     0,
		  0 ,     0, VAR_PZ,     0,     0,     0,
		  0 ,     0,     0, VAR_PYAW,     0,     0,
		  0 ,     0,     0,     0, VAR_PPITCH,     0,
		  0 ,     0,     0,     0,     0, VAR_PROLL;

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
		sigma(3, c) = sigma_x(9, c);			// yaw
		sigma(4, c) = sigma_x(10, c);			// pitch
		sigma(5, c) = sigma_x(11, c);			// roll
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
		dz(3) = normalize(dz(3));
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
