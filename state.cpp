#include "state.h"
#include "utils.h"
#include <iostream>

#include "m0.h"
#include "m1.h"
#include "m2.h"
#include "m3.h"

state_predict::state_predict() {
	_x.resize(NX);
	_x = VectorXd::Zero(NX);
	_P.resize(NX, NX);
	_P = MatrixXd::Identity(NX, NX);
}

void state_predict::initialize(int NSIGMA, int NAUG, double W, double W0_m, double W0_c, vector<double> noises, double SCALE) {

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

	_sigma.resize(NX, NSIGMA);

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

	/*
	BDCSVD<MatrixXd> USV(augmented_P);
	VectorXd diag = USV.singularValues();
	MatrixXd S = MatrixXd(diag.size(), diag.size());
	for (int i = 0; i < diag.size(); ++i) {
		if (diag(diag.size() - i - 1) < 10 ^ -10)
			S(i, i) = 1e-8;
		else
			S(i, i) = diag(diag.size() - i - 1);
		//S(i, i) = ((diag(diag.size() - i - 1) < 0) ? (1e-8) : (diag(diag.size() - i - 1)));
	}
	augmented_P = USV.computeU() * S * USV.computeV();
	*/

	MatrixXd L = augmented_P.llt().matrixL();
	//MatrixXd L = augmented_P.ldlt().matrixL();

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
	cerr << "UKF ERROR: DERIVED STATE FUNCTION FOR PREDICT_SIGMA IS NOT SET\n";
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
		predicted_P += _WEIGHTS_c[c] * dx * dx.transpose();
	}

	return predicted_P;
}

void state_predict::process(const VectorXd& current_x, const MatrixXd& current_P, const double dt) {

	MatrixXd augmented_sigma = compute_augmented_sigma(current_x, current_P);
	_sigma = predict_sigma(augmented_sigma, dt);
	_x = predict_x(_sigma);
	_P = predict_P(_sigma, _x);
}

VectorXd state_predict::peek(const VectorXd& current_x, const MatrixXd& current_P, double dt) {
	MatrixXd augmented_sigma = compute_augmented_sigma(current_x, current_P);
	MatrixXd sigma = predict_sigma(augmented_sigma, dt);
	return predict_x(sigma);
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
