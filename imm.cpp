#include "imm.h"
#include "m0.h"
#include "m1.h"
#include "m2.h"
#include "m3.h"

#include <iostream>

imm::imm(bool debug) {
	_x.resize(NX);
	_x = VectorXd::Zero(NX);
	_T.resize(NM, NM);
	_mu.resize(NM);
	_llh.resize(NM);
	_llh = VectorXd::Zero(NM);
	_normalizer = VectorXd::Zero(NM);
	_omega.resize(NM, NM);
	_X.resize(NX, NM);
	_X = MatrixXd::Zero(NX, NM);
	_p.resize(NX, NX);
	_P.resize(NM);
	for (int i = 0; i < NM; ++i)
		_P[i] = MatrixXd::Identity(NX, NX);

	_debug = debug;
}

imm::~imm() {
	for (int i = 0; i < NM; ++i) {
		delete _ukf[i];
	}
}

void imm::initialize(ALL_DATA, double t) {
	_T << 0.7, 0.1, 0.1, 0.1,
		0.1, 0.7, 0.1, 0.1,
		0.0, 1.0, 0.0, 0.0,
		0.1, 0.1, 0.1, 0.7;
	_mu << 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0;

	_ukf[0] = new ukf0(DATA, t, _debug);
	_ukf[1] = new ukf1(DATA, t, _debug);
	_ukf[2] = new ukf2(DATA, t, _debug);
	_ukf[3] = new ukf3(DATA, t, _debug);

	compute_mix_prob();

	_initialized = true;
}

void imm::mix(MatrixXd models_x) {

	// Mixing estimate and covariance
	MatrixXd new_x = MatrixXd::Zero(NX, NM);

	for (int i = 0; i < NM; ++i) {
		for (int j = 0; j < NM; ++j) {
			new_x.col(i) += models_x.col(j) * _omega(j, i);
		}
		_X.col(i) = new_x.col(i);
		MatrixXd new_P = MatrixXd::Zero(NX, NX);
		for (int j = 0; j < NM; ++j) {
			new_P += _omega(j, i) * (_P[j] + (models_x.col(j) - _X.col(i)) * (models_x.col(j) - _X.col(i)).transpose());
		}
		_P[i] = new_P;
	}
}

void imm::compute_state() {
	_x = VectorXd::Zero(NX);

	for (int i = 0; i < NM; ++i) {
		_x += (*_ukf[i]).get() * _mu(i);
	}

	_p = MatrixXd::Zero(NX, NX);

	for (int i = 0; i < NM; ++i) {
		_p += _mu(i) * (_P[i] + ((*_ukf[i]).get() - _x) * ((*_ukf[i]).get() - _x).transpose());
	}
}

void imm::compute_mix_prob() {
	_normalizer = _mu.adjoint() * _T;

	// Predicted model probability
	for (int i = 0; i < NM; ++i) {
		for (int j = 0; j < NM; ++j) {
			_omega(i, j) = (_T(i, j) * _mu(i)) / _normalizer(j);
		}
	}
}

void imm::update(ALL_DATA, double t) {

	for (int i = 0; i < NM; ++i) {
		(*_ukf[i]).update(DATA, _X.col(i), _P[i]);
		_llh(i) = (*_ukf[i]).get_llh();
		//cout << _llh(i) << " ";
	}
	
	double mu_sum = 0;

	for (int i = 0; i < NM; ++i) {
		_mu(i) = _normalizer(i) * _llh(i);
		mu_sum += _mu(i);
	}

	for (int i = 0; i < NM; ++i) {
		_mu(i) /= mu_sum;
	}

	compute_mix_prob();
	compute_state();
}

void imm::compute(ALL_DATA, double t) {
	
	predict(t);
	update(DATA, t);
}

void imm::process(ALL_DATA, double t) {
	_initialized ? compute(DATA, t) : initialize(DATA, t);
}

void imm::predict(double t) {

	MatrixXd models_x(NX, NM);
	for (int i = 0; i < NM; ++i) {
		models_x.col(i) = (*_ukf[i]).get();
	}

	mix(models_x);

	for (int i = 0; i < NM; ++i) {
		(*_ukf[i]).predict(t, _X.col(i), _P[i]);
	}

	compute_state();
}

VectorXd imm::peek(double dt) {

	MatrixXd ukf_x(NX, NM);
	VectorXd peek_x = VectorXd::Zero(NX);

	for (int i = 0; i < NM; ++i) {
		ukf_x.col(i) = (*_ukf[i]).peek(dt);
	}

	for (int i = 0; i < NM; ++i) {
		peek_x += ukf_x.col(i) * _mu(i);
	}

	return peek_x;
}

VectorXd imm::get() const {
	return _x;
}

VectorXd imm::get_mu() {
	return _mu;
}

VectorXd imm::get_xi(int i) {
	return (*_ukf[i]).get();
}

VectorXd imm::get_xi_peek(int i, double dt) {
	return (*_ukf[i]).peek(dt);
}
