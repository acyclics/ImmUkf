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

	_R(3, 3) = VAR_PYAW;
	_R(4, 4) = VAR_PPITCH;
	_R(5, 5) = VAR_PROLL;

	_R(6, 6) = VAR_XVEL;
	_R(7, 7) = VAR_YVEL;
	_R(8, 8) = VAR_ZVEL;

	_R(9, 9) = VAR_YAWVEL;
	_R(10, 10) = VAR_PITCHVEL;
	_R(11, 11) = VAR_ROLLVEL;

	_R(12, 12) = VAR_XACC;
	_R(13, 13) = VAR_YACC;
	_R(14, 14) = VAR_ZACC;

	_R(15, 15) = VAR_YAWACC;
	_R(16, 16) = VAR_PITCHACC;
	_R(17, 17) = VAR_ROLLACC;

	_R(18, 18) = VAR_DIST2CENTER;

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

		sigma(6, c) = sigma_x(3, c);			// xvel
		sigma(7, c) = sigma_x(4, c);			// yvel
		sigma(8, c) = sigma_x(5, c);			// zvel

		sigma(9, c) = sigma_x(12, c);			// yawvel
		sigma(10, c) = sigma_x(13, c);			// pitchvel
		sigma(11, c) = sigma_x(14, c);			// rollvel

		sigma(12, c) = sigma_x(6, c);			// xacc
		sigma(13, c) = sigma_x(7, c);			// yacc
		sigma(14, c) = sigma_x(8, c);			// zacc

		sigma(15, c) = sigma_x(15, c);			// yawacc
		sigma(16, c) = sigma_x(16, c);			// pitchacc
		sigma(17, c) = sigma_x(17, c);			// rollacc

		sigma(18, c) = sigma_x(18, c);			// dist2center
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
