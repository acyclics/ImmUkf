#include "imm.h"
#include "m0.h"
#include "m1.h"

imm::imm() {
}

void imm::initialize(ALL_DATA, double t) {
	_T << 0.7, 0.3,
		  0.5, 0.5;
	_u << 0.5, 0,
		  0.5, 0;

	_ukf[0] = ukf0(px, py, pz, yaw, pitch, roll, t);
	_ukf[1] = ukf1(px, py, pz, yaw, pitch, roll, t);
}

void imm::mix(MatrixXd models_x) {

	MatrixXd new_u = MatrixXd::Zero(NM, NM);
	double normalizer = 0;

	// Predicted model probability
	for (int i = 0; i < NM; ++i) {
		for (int j = 0; j < NM; ++j) {
			new_u(i, i) += _T(j, i) * _u(j);
		}
		_u(i, i) = new_u(i, i);
		normalizer += _u(i, i);
	}

	// Mixing probability
	for (int i = 0; i < NM; ++i) {
		for (int j = 0; j < NM; ++j) {
			_u(j, i) = _T(j, i) * _u(j, j) / normalizer;
		}
	}

	// Mixing estimate
	MatrixXd new_x = MatrixXd::Zero(NX, NM);

	for (int i = 0; i < NM; ++i) {
		for (int j = 0; j < NM; ++j) {
			new_x.col(i) += models_x.col(j) * _u(j, i);
		}
		_X.col(i) = new_x.col(i);
	}

	// Mixing covariance
	for (int i = 0; i < NM; ++i) {
		MatrixXd new_P = MatrixXd::Zero(NM, NM);
		for (int j = 0; j < NM; ++j) {
			new_P += _u(j, i) * (_P[j] + (models_x.col(j) - _X.col(i)) * (models_x.col(j) - _X.col(i)).transpose());
		}
		_P[i] = new_P;
	}
}

void imm::filtering(ALL_DATA, double t) {
	for (int i = 0; i < NM; ++i) {
		_ukf[i].process(px, py, pz, yaw, pitch, roll, t, _X.col(i), _P[i]);
	}
}

void imm::compute(ALL_DATA, double t) {

	MatrixXd models_x(NX, NM);
	for (int i = 0; i < NM; ++i) {
		models_x.col(i) = _ukf[i].get();
	}

	mix(models_x);
	filtering(px, py, pz, yaw, pitch, roll, t);

	VectorXd new_x = VectorXd::Zero(NX);
	for (int i = 0; i < NM; ++i) {
		new_x += _ukf[i].get() * _u(i, i);
	}
	_x = new_x;
}

void imm::process(ALL_DATA, double t) {
	_initialized ? compute(px, py, pz, yaw, pitch, roll, t) : initialize(px, py, pz, yaw, pitch, roll, t);
}
