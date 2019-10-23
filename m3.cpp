#include "m3.h"
#include <iostream>
#include "global_settings.h"
#include "ukf.h"
#include "m3_settings.h"
#include "utils.h"

model3::model3() {
}

void model3::update(double vx, double vy, double vz, double qw, double qx, double qy, double qz, double new_measurement_t) {

	Quaterniond q(qw, qx, qy, qz);
	q = q.normalized();
	auto rollpitchyaw = q.toRotationMatrix().eulerAngles(0, 1, 2);

	_current_vxyz(0) = vx;
	_current_vxyz(1) = vy;
	_current_vxyz(2) = vz;

	_measure_rollpitchyaw(0) = rollpitchyaw[0];
	_measure_rollpitchyaw(1) = rollpitchyaw[1];
	_measure_rollpitchyaw(2) = rollpitchyaw[2];

	_prev_measurement_t = _measurement_t;
	_measurement_t = new_measurement_t;
}

Vector3d model3::get_current_vxyz() {
	return _current_vxyz;
}

Vector3d model3::get_measure_rollpitchyaw() {
	return _measure_rollpitchyaw;
}

state3::state3() : state_predict() {
}

MatrixXd state3::predict_sigma(const MatrixXd& augmented_sigma, double dt) {

	MatrixXd predicted_sigma = MatrixXd(NX, _NSIGMA);

	for (int c = 0; c < _NSIGMA; ++c) {
		/*
			1. Get the current state
		*/
		const double px = augmented_sigma(0, c);
		const double py = augmented_sigma(1, c);
		double pz = augmented_sigma(2, c);

		const double proll = augmented_sigma(3, c);
		const double ppitch = augmented_sigma(4, c);
		const double pyaw = augmented_sigma(5, c);

		const double vx_noise = augmented_sigma(6, c);
		const double vy_noise = augmented_sigma(7, c);
		const double vz_noise = augmented_sigma(8, c);

		const double proll_noise = augmented_sigma(9, c);
		const double ppitch_noise = augmented_sigma(10, c);
		const double pyaw_noise = augmented_sigma(11, c);

		/*
			2. Predict next state with noise
		*/
		double dt2 = dt * dt;
		
		predicted_sigma(0, c) = px + ((model3*)_m)->get_current_vxyz()(0) * dt + vx_noise;
		predicted_sigma(1, c) = py + ((model3*)_m)->get_current_vxyz()(1) * dt + vy_noise;
		predicted_sigma(2, c) = pz + ((model3*)_m)->get_current_vxyz()(2) * dt + vz_noise;

		predicted_sigma(3, c) = fmodl(((model3*)_m)->get_measure_rollpitchyaw()(0) + proll_noise, 2 * M_PI);
		predicted_sigma(4, c) = fmodl(((model3*)_m)->get_measure_rollpitchyaw()(1) + ppitch_noise, 2 * M_PI);
		predicted_sigma(5, c) = fmodl(((model3*)_m)->get_measure_rollpitchyaw()(2) + pyaw_noise, 2 * M_PI);
	}

	return predicted_sigma;
}

ukf3::ukf3(ALL_DATA, double t, model* m, bool debug) : ukf() {
	using namespace m3_sets;
	using namespace std;
	vector<double> noises{ VAR_VX_NOISE, VAR_VY_NOISE, VAR_VZ_NOISE, VAR_PROLL_NOISE, VAR_PPITCH_NOISE, VAR_PYAW_NOISE };
	_P = MatrixXd::Identity(NX, NX) * 1.0;

	initialize(DATA, t, NSIGMA, NAUG, W, W0_m, W0_c, noises, SCALE, VAR, m, debug);
}

state_predict& ukf3::get_stater() {
	return _stater3;
}
