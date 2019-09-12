#include "m3.h"
#include <iostream>
#include "global_settings.h"
#include "ukf.h"
#include "m3_settings.h"

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
		const double pz = augmented_sigma(2, c);

		const double vx = augmented_sigma(3, c);
		const double vy = augmented_sigma(4, c);
		const double vz = augmented_sigma(5, c);

		const double ax = augmented_sigma(6, c);
		const double ay = augmented_sigma(7, c);
		const double az = augmented_sigma(8, c);

		const double yaw = augmented_sigma(9, c);
		const double pitch = augmented_sigma(10, c);
		const double roll = augmented_sigma(11, c);

		const double vyaw = augmented_sigma(12, c);
		const double vpitch = augmented_sigma(13, c);
		const double vroll = augmented_sigma(14, c);

		const double ayaw = augmented_sigma(15, c);
		const double apitch = augmented_sigma(16, c);
		const double aroll = augmented_sigma(17, c);

		const double dist2center = augmented_sigma(18, c);

		const double vx_noise = augmented_sigma(19, c);
		const double vy_noise = augmented_sigma(20, c);
		const double vz_noise = augmented_sigma(21, c);

		const double pyaw_noise = augmented_sigma(22, c);
		const double ppitch_noise = augmented_sigma(23, c);
		const double proll_noise = augmented_sigma(24, c);

		/*
			2. Predict next state with noise
		*/
		double dt2 = dt * dt;

		predicted_sigma(0, c) = px + vx * dt;
		predicted_sigma(1, c) = py + vy * dt;
		predicted_sigma(2, c) = pz + vz * dt;

		predicted_sigma(3, c) = vx + vx_noise;
		predicted_sigma(4, c) = vy + vy_noise;
		predicted_sigma(5, c) = vz + vz_noise;

		predicted_sigma(6, c) = 0;
		predicted_sigma(7, c) = 0;
		predicted_sigma(8, c) = 0;

		predicted_sigma(9, c) = yaw + pyaw_noise;
		predicted_sigma(10, c) = pitch + ppitch_noise;
		predicted_sigma(11, c) = roll + proll_noise;

		predicted_sigma(12, c) = 0;// vyaw;
		predicted_sigma(13, c) = 0;// vpitch;
		predicted_sigma(14, c) = 0;// vroll;

		predicted_sigma(15, c) = 0;// ayaw;
		predicted_sigma(16, c) = 0;// apitch;
		predicted_sigma(17, c) = 0;// aroll;

		predicted_sigma(18, c) = dist2center;
	}

	return predicted_sigma;
}

ukf3::ukf3(ALL_DATA, double t, bool debug) : ukf() {
	using namespace m3_sets;
	using namespace std;
	vector<double> noises{ VAR_VX_NOISE, VAR_VY_NOISE, VAR_VZ_NOISE, VAR_PYAW_NOISE, VAR_PPITCH_NOISE, VAR_PROLL_NOISE };
	_P = MatrixXd::Identity(NX, NX) * 1.0;

	/*
	_P(6, 6) = 0;
	_P(7, 7) = 0;
	_P(8, 8) = 0;

	_P(12, 12) = 0;
	_P(13, 13) = 0;
	_P(14, 14) = 0;
	*/

	_P(6, 6) = 0;
	_P(7, 7) = 0;
	_P(8, 8) = 0;

	_P(12, 12) = 0;
	_P(13, 13) = 0;
	_P(14, 14) = 0;

	_P(15, 15) = 0;
	_P(16, 16) = 0;
	_P(17, 17) = 0;

	_P(18, 18) = 0;

	initialize(DATA, t, NSIGMA, NAUG, W, W0_m, W0_c, noises, SCALE, VAR, debug);
}

state_predict& ukf3::get_stater() {
	return _stater3;
}
