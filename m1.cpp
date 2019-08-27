#include "m1.h"
#include "m1_settings.h"

ukf ukf1(ALL_DATA, double t) {
	using namespace m1_sets;
	using namespace std;
	ukf ukf1;
	vector<double> noises = { VAR_VX_NOISE, VAR_VY_NOISE, VAR_VZ_NOISE, VAR_VYAW_NOISE, VAR_VPITCH_NOISE, VAR_VROLL_NOISE, VAR_DIST2CENTER_NOISE };
	ukf1.initialize(px, py, pz, yaw, pitch, roll, t, 0, NSIGMA, NAUG, W, W0_m, W0_c, noises, SCALE, VAR_PX, VAR_PY, VAR_PZ, VAR_PYAW, VAR_PPITCH, VAR_PROLL);
	return ukf1;
}

void model1(MatrixXd& predicted_sigma, const MatrixXd& augmented_sigma, int NSIGMA, double dt) {
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

		const double yaw = augmented_sigma(9, c);
		const double pitch = augmented_sigma(10, c);
		const double roll = augmented_sigma(11, c);

		const double vyaw = augmented_sigma(12, c);
		const double vpitch = augmented_sigma(13, c);
		const double vroll = augmented_sigma(14, c);

		const double dist2center = augmented_sigma(18, c);

		const double vx_noise = augmented_sigma(19, c);
		const double vy_noise = augmented_sigma(20, c);
		const double vz_noise = augmented_sigma(21, c);

		const double vyaw_noise = augmented_sigma(22, c);
		const double vpitch_noise = augmented_sigma(23, c);
		const double vroll_noise = augmented_sigma(24, c);

		const double dist2center_noise = augmented_sigma(25, c);

		/*
			2. Predict next state with noise
		*/
		predicted_sigma(0, c) = px + vx * dt + dist2center * vyaw * cos(yaw) * dt + dist2center * vpitch * cos(pitch) * dt;
		predicted_sigma(1, c) = py + vy * dt + dist2center * vyaw * sin(yaw) * dt + dist2center * vroll * cos(roll) * dt;
		predicted_sigma(2, c) = pz + vz * dt  + dist2center * vroll * sin(roll) * dt + dist2center * vpitch * sin(pitch) * dt;

		predicted_sigma(3, c) = vx + vx_noise * dt;
		predicted_sigma(4, c) = vy + vy_noise * dt;
		predicted_sigma(5, c) = vz + vz_noise * dt;

		predicted_sigma(6, c) = 0;
		predicted_sigma(7, c) = 0;
		predicted_sigma(8, c) = 0;

		predicted_sigma(9, c) = yaw + vyaw * dt;
		predicted_sigma(10, c) = pitch + vpitch * dt;
		predicted_sigma(11, c) = roll + vroll * dt;

		predicted_sigma(12, c) = vyaw + vyaw_noise * dt;
		predicted_sigma(13, c) = vpitch + vpitch_noise * dt;
		predicted_sigma(14, c) = vroll + vroll_noise * dt;

		predicted_sigma(15, c) = 0;
		predicted_sigma(16, c) = 0;
		predicted_sigma(17, c) = 0;

		predicted_sigma(18, c) = dist2center + dist2center_noise;
	}
}
