#include "ukf.h"

using namespace Eigen;

ukf::ukf() {
	_initialized = false;
}

void ukf::initialize(ALL_DATA, double t, int model_no, int NSIGMA, int NAUG, double W, double W0_m, double W0_c, vector<double> noises, double SCALE,
					 double VAR_PX, double VAR_PY, double VAR_PZ, double VAR_PYAW, double VAR_PPITCH, double VAR_PROLL) {
	_x = VectorXd(6);
	_x << px, py, pz, yaw, pitch, roll;
	_P = MatrixXd::Identity(NX, NX);
	_t = t;
	_initialized = true;

	_model_no = model_no;

	_stater.initialize(model_no, NSIGMA, NAUG, W, W0_m, W0_c, noises, SCALE);
	_measurer.initialize(NSIGMA, W, W0_m, W0_c, VAR_PX, VAR_PY, VAR_PZ, VAR_PYAW, VAR_PPITCH, VAR_PROLL);
	_merger.initialize(NSIGMA, W, W0_m, W0_c);
}

void ukf::update(ALL_DATA, double t, VectorXd imm_x, MatrixXd imm_P) {

	VectorXd predicted_z;
	MatrixXd sigma_x;
	MatrixXd sigma_z;
	MatrixXd S;

	VectorXd z = VectorXd(6);
	z << px, py, pz, yaw, pitch, roll;

	double dt = (t - _t);

	// STATE PREDICTION
	_stater.process(imm_x, imm_P, dt);
	_x = _stater.getx();
	_P = _stater.getP();
	sigma_x = _stater.get_sigma();

	// MEASUREMENT PREDICTION
	_measurer.process(sigma_x);
	predicted_z = _measurer.getz();
	S = _measurer.getS();
	sigma_z = _measurer.get_sigma();

	// STATE UPDATE
	_merger.process(_x, predicted_z, z, S, _P, sigma_x, sigma_z);
	_x = _merger.getx();
	_P = _merger.getP();
	_nis = _merger.get_nis();

	// update t
	_t = t;
}

void ukf::process(ALL_DATA, double t, VectorXd imm_x, MatrixXd imm_P) {
	update(px, py, pz, yaw, pitch, roll, t, imm_x, imm_P);
}

VectorXd ukf::get() const {
	return _x;
}

double ukf::get_nis() const {
	return _nis;
}
