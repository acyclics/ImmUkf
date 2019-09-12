#include "m2.h"
#include "m2_settings.h"
#include "utils.h"
#include <iostream>

model2 m2;

state2::state2() : state_predict() {
}

MatrixXd state2::predict_sigma(const MatrixXd& augmented_sigma, double dt) {

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

		const double px_noise = augmented_sigma(19, c);
		const double py_noise = augmented_sigma(20, c);
		const double pz_noise = augmented_sigma(21, c);

		const double pyaw_noise = augmented_sigma(22, c);
		const double ppitch_noise = augmented_sigma(23, c);
		const double proll_noise = augmented_sigma(24, c);

		const double vx_noise = augmented_sigma(25, c);
		const double vy_noise = augmented_sigma(26, c);
		const double vz_noise = augmented_sigma(27, c);

		const double vyaw_noise = augmented_sigma(28, c);
		const double vpitch_noise = augmented_sigma(29, c);
		const double vroll_noise = augmented_sigma(30, c);

		const double ax_noise = augmented_sigma(31, c);
		const double ay_noise = augmented_sigma(32, c);
		const double az_noise = augmented_sigma(33, c);

		const double dist2center_noise = augmented_sigma(34, c);

		/*
			2. Predict next state with noise
		*/

		Quaterniond old_q = AngleAxisd(roll, Vector3d::UnitX()) * AngleAxisd(pitch, Vector3d::UnitY()) * AngleAxisd(yaw, Vector3d::UnitZ());
		old_q = old_q.normalized();
		Vector3d w(vroll, vpitch, vyaw);
		Vector3d pos(px, py, pz);

		Eigen::MatrixXd J_e = MatrixXd(3, 4);
		J_e << -m2.q.x(), m2.q.w(), -m2.q.z(), m2.q.y(), -m2.q.y(), m2.q.z(), m2.q.w(), -m2.q.x(), -m2.q.z(), -m2.q.y(), m2.q.x(), m2.q.w();
		Vector4d qdiff = Eigen::Vector4d(m2.q.w() - m2.prev_q.w(), m2.q.x() - m2.prev_q.x(), m2.q.y() - m2.prev_q.y(), m2.q.z() - m2.prev_q.z());
		Vector3d vel_rollpitchyaw = vel_rollpitchyaw = 2 * J_e * qdiff;

		Vector3d u = pos - m2.h;

		Vector3d vxyz = w.cross(u);
		Vector3d axyz = w.cross(w.cross(u));

		Quaterniond d_q = (old_q * Quaterniond(0, vroll, vpitch, vyaw));
		d_q.coeffs() *= 0.5;

		Quaterniond new_q;
		new_q.coeffs() = old_q.coeffs() + d_q.coeffs() * dt;
		new_q = new_q.normalized();
		auto new_orien = new_q.toRotationMatrix().eulerAngles(0, 1, 2);

		Vector3d new_pos = (AngleAxisd(w.norm() * dt, w.normalized()).toRotationMatrix() * u) + m2.h;

		predicted_sigma(0, c) = new_pos[0] + px_noise;
		predicted_sigma(1, c) = new_pos[1] + py_noise;
		predicted_sigma(2, c) = new_pos[2] + pz_noise;

		predicted_sigma(3, c) = vxyz[0] + vx_noise;
		predicted_sigma(4, c) = vxyz[1] + vy_noise;
		predicted_sigma(5, c) = vxyz[2] + vz_noise;

		predicted_sigma(6, c) = axyz[0] + ax_noise;
		predicted_sigma(7, c) = axyz[1] + ay_noise;
		predicted_sigma(8, c) = axyz[2] + az_noise;

		predicted_sigma(9, c) = new_orien[2] + pyaw_noise;
		predicted_sigma(10, c) = new_orien[1] + ppitch_noise;
		predicted_sigma(11, c) = new_orien[0] + proll_noise;

		predicted_sigma(12, c) = vel_rollpitchyaw[2] + vyaw_noise;
		predicted_sigma(13, c) = vel_rollpitchyaw[1] + vpitch_noise;
		predicted_sigma(14, c) = vel_rollpitchyaw[0] + vroll_noise;

		predicted_sigma(15, c) = ayaw;
		predicted_sigma(16, c) = apitch;
		predicted_sigma(17, c) = aroll;

		predicted_sigma(18, c) = dist2center + dist2center_noise;
	}

	return predicted_sigma;
}

ukf2::ukf2(ALL_DATA, double t, bool debug) : ukf() {
	using namespace m2_sets;
	using namespace std;
	vector<double> noises{ VAR_PX_NOISE, VAR_PY_NOISE, VAR_PZ_NOISE, VAR_PYAW_NOISE, VAR_PPITCH_NOISE, VAR_PROLL_NOISE, VAR_VX_NOISE, VAR_VY_NOISE, VAR_VZ_NOISE, VAR_VYAW_NOISE, VAR_VPITCH_NOISE, VAR_VROLL_NOISE, VAR_AX_NOISE, VAR_AY_NOISE, VAR_AZ_NOISE, VAR_DIST2CENTER_NOISE };
	_P = MatrixXd::Identity(NX, NX) * 10;
	/*
	_P(0, 0) = 0;
	_P(1, 1) = 0;
	_P(2, 2) = 0;
	_P(3, 3) = 0;
	_P(4, 4) = 0;
	_P(5, 5) = 0;
	_P(6, 6) = 0;
	_P(7, 7) = 0;
	_P(8, 8) = 0;*/

	
	_P(15, 15) = 0;
	_P(16, 16) = 0;
	_P(17, 17) = 0;
	
	initialize(DATA, t, NSIGMA, NAUG, W, W0_m, W0_c, noises, SCALE, VAR, debug);
}

state_predict& ukf2::get_stater() {
	return _stater2;
}
