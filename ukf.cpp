#include "ukf.h"
#include "utils.h"
#include <iostream>

using namespace Eigen;

ukf::ukf() {
	_initialized = false;
	_x = VectorXd::Zero(NX);
}

void ukf::initialize(ALL_DATA, double t, int NSIGMA, int NAUG, double W, double W0_m, double W0_c, vector<double> noises, double SCALE,
					 ALL_VAR, bool debug) {
	
	_t = t;
	_initialized = true;
	_debug = debug;
	_noises = noises;

	get_stater().initialize(NSIGMA, NAUG, W, W0_m, W0_c, noises, SCALE);
	_measurer.initialize(NSIGMA, W, W0_m, W0_c, VAR);
	_merger.initialize(NSIGMA, W, W0_m, W0_c);
}

void ukf::process(ALL_DATA, double t, const VectorXd& imm_x, const MatrixXd& imm_P) {
	
	predict(t, imm_x, imm_P);
	update(DATA, imm_x, imm_P);
}

void ukf::predict(double t, const VectorXd& imm_x, const MatrixXd& imm_P) {

	double dt = (t - _t);

	// STATE PREDICTION
	if (_debug)
		get_stater().process(_x, _P, dt);
	else
		get_stater().process(imm_x, imm_P, dt);
	_x = get_stater().getx();
	_P = get_stater().getP();
	//_P = pdefinite_svd(_P);
	_sigma_x = get_stater().get_sigma();

	// MEASUREMENT PREDICTION
	_measurer.process(_sigma_x);
	_predicted_z = _measurer.getz();
	_S = _measurer.getS();
	//_S = pdefinite_svd(_S);
	_sigma_z = _measurer.get_sigma();

	// update t
	_t = t;
}

void ukf::update(ALL_DATA, const VectorXd& imm_x, const MatrixXd& imm_P) {

	VectorXd z = VectorXd(NZ);
	z << DATA;

	// STATE UPDATE
	if (_debug)
		_merger.process(_x, _predicted_z, z, _S, _P, _sigma_x, _sigma_z);
	else
		_merger.process(_x, _predicted_z, z, _S, _P, _sigma_x, _sigma_z);
		//_merger.process(imm_x, _predicted_z, z, _S, imm_P, _sigma_x, _sigma_z);
	_x = _merger.getx();
	_P = _merger.getP();
	_nis = _merger.get_nis();

	// calculate likelihood
	//_llh = (1 / (sqrtl(2 * M_PI * _S.determinant()))) * (expl((-1.0 / 2.0) * _predicted_z.transpose() * _S.inverse() * _predicted_z));
	double log_llh = logl(1) - 0.5 * (logl(2) + logl(M_PI) + logl(_S.determinant())) + (-0.5 * _predicted_z.transpose() * _S.inverse() * _predicted_z);
	//cout << "Printing z\n" << _S << "\n";
	_llh = expl(log_llh);
	//cout << "Printing llh\n" << log_llh << "\n";
	//_llh = logl(labs(_tmp_llh)) * labs(_tmp_llh) / (_tmp_llh);
}

VectorXd ukf::peek(double dt) {
	VectorXd peek_x = get_stater().peek(_x, _P, dt);
	return peek_x;
}

VectorXd ukf::get() const {
	return _x;
}

MatrixXd ukf::get_P() const {
	return _P;
}

double ukf::get_nis() const {
	return _nis;
}

double ukf::get_llh() const {
	return _llh;
}

double ukf::get_measurement_pred(int no) const {
	return _predicted_z(no);
}

double ukf::get_state_var(int i, int j) {
	return _P(i, j);
}

double ukf::get_noise(int i) const {
	return _noises[i];
}

state_predict& ukf::get_stater() {
	state_predict placeholder;
	cerr << "UKF ERROR: DERIVED UKF FUNCTION FOR GET_STATER IS NOT SET\n";
	return placeholder;
}
