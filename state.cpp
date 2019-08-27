#include "state.h"
#include "utils.h"

#include "m0.h"
#include "m1.h"

state_predict::state_predict() {
}

void state_predict::initialize(int model_no, int NSIGMA, int NAUG, double W, double W0_m, double W0_c, vector<double> noises, double SCALE) {

	_model_no = model_no;
	_NSIGMA = NSIGMA;
	_NAUG = NAUG;

	_WEIGHTS_m.resize(NSIGMA);
	_WEIGHTS_m[0] = W0_m;
	_WEIGHTS_c.resize(NSIGMA);
	_WEIGHTS_c[0] = W0_c;

	for (int i = 1; i < NSIGMA; ++i) {
		_WEIGHTS_m[i] = W;
		_WEIGHTS_c[i] = W;
	}

	_sigma = MatrixXd(NX, NSIGMA);

	_noises = noises;
	_noises_size = noises.size();

	_SCALE = SCALE;
}

MatrixXd state_predict::compute_augmented_sigma(const VectorXd& current_x, const MatrixXd& current_P) {

	MatrixXd augmented_sigma = MatrixXd::Zero(_NAUG, _NSIGMA);
	VectorXd augmented_x = VectorXd::Zero(_NAUG);
	MatrixXd augmented_P = MatrixXd::Zero(_NAUG, _NAUG);

	augmented_x.head(NX) = current_x;

	augmented_P.topLeftCorner(NX, NX) = current_P;
	for (int i = 0; i < _noises_size; ++i) {
		augmented_P(NX + i, NX + i) = _noises[i];
	}

	const MatrixXd L = augmented_P.llt().matrixL();
	augmented_sigma.col(0) = augmented_x;

	for (int c = 0; c < _NAUG; ++c) {
		const int i = c + 1;
		augmented_sigma.col(i) = augmented_x + _SCALE * L.col(c);
		augmented_sigma.col(i + _NAUG) = augmented_x - _SCALE * L.col(c);
	}

	return augmented_sigma;
}

MatrixXd state_predict::predict_sigma(const MatrixXd& augmented_sigma, double dt) {

	MatrixXd predicted_sigma = MatrixXd(NX, _NSIGMA);

	switch (_model_no) {
		case 0:
			model0(predicted_sigma, augmented_sigma, _NSIGMA, dt);
			break;
		
		case 1:
			model1(predicted_sigma, augmented_sigma, _NSIGMA, dt);
			break;
	}

	return predicted_sigma;
}

VectorXd state_predict::predict_x(const MatrixXd& predicted_sigma) {

	VectorXd predicted_x = VectorXd::Zero(NX);

	for (int c = 0; c < _NSIGMA; ++c) {
		predicted_x += _WEIGHTS_m[c] * predicted_sigma.col(c);
	}

	return predicted_x;
}

MatrixXd state_predict::predict_P(const MatrixXd& predicted_sigma, const VectorXd& predicted_x) {

	MatrixXd predicted_P = MatrixXd::Zero(NX, NX);
	VectorXd dx = VectorXd(NX);

	for (int c = 0; c < _NSIGMA; c++) {

		dx = predicted_sigma.col(c) - predicted_x;
		dx(3) = normalize(dx(3));
		predicted_P += _WEIGHTS_c[c] * dx * dx.transpose();
	}

	return predicted_P;
}

void state_predict::process(VectorXd& current_x, MatrixXd& current_P, double dt) {

	MatrixXd augmented_sigma = compute_augmented_sigma(current_x, current_P);
	_sigma = predict_sigma(augmented_sigma, dt);
	_x = predict_x(_sigma);
	_P = predict_P(_sigma, _x);
}

MatrixXd state_predict::getP() const {
	return _P;
}

MatrixXd state_predict::get_sigma() const {
	return _sigma;
}

VectorXd state_predict::getx() const {
	return _x;
}
