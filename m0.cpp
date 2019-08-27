#include "m0.h"
#include "m0_settings.h"

ukf ukf0(ALL_DATA, double t) {
	using namespace m0_sets;
	using namespace std;
	ukf ukf0;
	vector<double> noises = { VAR_AX_NOISE, VAR_AY_NOISE, VAR_AZ_NOISE, VAR_AYAW_NOISE, VAR_APITCH_NOISE, VAR_AROLL_NOISE, VAR_DIST2CENTER_NOISE };
	ukf0.initialize(px, py, pz, yaw, pitch, roll, t, 0, NSIGMA, NAUG, W, W0_m, W0_c, noises, SCALE, VAR_PX, VAR_PY, VAR_PZ, VAR_PYAW, VAR_PPITCH, VAR_PROLL);
	return ukf0;
}

void model0(MatrixXd& predicted_sigma, const MatrixXd& augmented_sigma, int NSIGMA, double dt) {
	for (int c = 0; c < NSIGMA; ++c) {
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

		const double ax_noise = augmented_sigma(19, c);
		const double ay_noise = augmented_sigma(20, c);
		const double az_noise = augmented_sigma(21, c);

		const double ayaw_noise = augmented_sigma(22, c);
		const double apitch_noise = augmented_sigma(23, c);
		const double aroll_noise = augmented_sigma(24, c);

		const double dist2center_noise = augmented_sigma(25, c);

		/*
			2. Predict next state with noise
		*/
		double dt2 = dt * dt;

		predicted_sigma(0, c) = px + vx * dt + 0.5 * ax * dt2 + dist2center * vyaw * cos(yaw) * dt + 0.5 * ayaw * cos(yaw) * dt2 + dist2center * vpitch * cos(pitch) * dt + 0.5 * apitch * cos(pitch) * dt2;
		predicted_sigma(1, c) = py + vy * dt + 0.5 * ay * dt2 + dist2center * vyaw * sin(yaw) * dt + 0.5 * ayaw * sin(yaw) * dt2 + dist2center * vroll * cos(roll) * dt + 0.5 * aroll * cos(roll) * dt2;
		predicted_sigma(2, c) = pz + vz * dt + 0.5 * az * dt2 + dist2center * vroll * sin(roll) * dt + 0.5 * aroll * sin(roll) * dt2 + dist2center * vpitch * sin(pitch) * dt + 0.5 * apitch * sin(pitch) * dt2;

		predicted_sigma(3, c) = vx + ax * dt;
		predicted_sigma(4, c) = vy + ay * dt;
		predicted_sigma(5, c) = vz + az * dt;

		predicted_sigma(6, c) = ax + ax_noise * dt;
		predicted_sigma(7, c) = ay + ay_noise * dt;
		predicted_sigma(8, c) = az + az_noise * dt;

		predicted_sigma(9, c) = yaw + vyaw * dt + 0.5 * ayaw * dt * dt;
		predicted_sigma(10, c) = pitch + vpitch * dt + 0.5 * apitch * dt * dt;
		predicted_sigma(11, c) = roll + vroll * dt + 0.5 * aroll * dt * dt;

		predicted_sigma(12, c) = vyaw + ayaw * dt;
		predicted_sigma(13, c) = vpitch + apitch * dt;
		predicted_sigma(14, c) = vroll + aroll * dt;

		predicted_sigma(15, c) = ayaw + ayaw_noise * dt;
		predicted_sigma(16, c) = apitch + apitch_noise * dt;
		predicted_sigma(17, c) = aroll + aroll_noise * dt;

		predicted_sigma(18, c) = dist2center + dist2center_noise;
	}
}
