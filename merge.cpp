#include "merge.h"
#include "utils.h"

merger::merger() {
}

void merger::initialize(int NSIGMA, double W, double W0_m, double W0_c) {

	_NSIGMA = NSIGMA;

	_WEIGHTS_m.resize(NSIGMA);
	_WEIGHTS_m[0] = W0_m;
	_WEIGHTS_c.resize(NSIGMA);
	_WEIGHTS_c[0] = W0_c;

	for (int i = 1; i < NSIGMA; ++i) {
		_WEIGHTS_m[i] = W;
		_WEIGHTS_c[i] = W;
	}
}

MatrixXd merger::compute_Tc(const VectorXd& predicted_x, const VectorXd& predicted_z, const MatrixXd& sigma_x, const MatrixXd& sigma_z) {

	VectorXd dz;
	VectorXd dx;
	MatrixXd Tc = MatrixXd::Zero(NX, NZ);

	for (int c = 0; c < _NSIGMA; c++) {
		dx = sigma_x.col(c) - predicted_x;
		dz = sigma_z.col(c) - predicted_z;
		Tc += _WEIGHTS_c[c] * dx * dz.transpose();
	}

	return Tc;
}

void merger::update(const VectorXd& z, const MatrixXd& S, const MatrixXd& Tc, const VectorXd& predicted_z,
					const VectorXd& predicted_x, const MatrixXd& predicted_P) {

	MatrixXd Si = S.inverse();
	MatrixXd K = Tc * Si;

	VectorXd dz = z - predicted_z;

	_x = predicted_x + K * dz;
	_P = predicted_P - K * S * K.transpose();
	_nis = dz.transpose() * Si * dz;
}

void merger::process(const VectorXd& predicted_x, const VectorXd& predicted_z, const VectorXd& z, const MatrixXd& S, const MatrixXd& predicted_P,
					 const MatrixXd& sigma_x, const MatrixXd& sigma_z) {

	MatrixXd Tc = compute_Tc(predicted_x, predicted_z, sigma_x, sigma_z);
	update(z, S, Tc, predicted_z, predicted_x, predicted_P);
}

VectorXd merger::getx() const {

	return _x;
}

MatrixXd merger::getP() const {

	return _P;
}

double merger::get_nis() const {

	return _nis;
}
