#include "m1.h"
#include "m1_settings.h"
#include "utils.h"
#include <iostream>

model1::model1() {
}

void model1::update(VectorXd& new_h, Quaterniond& new_q, double new_measurement_t) {

	_h = new_h;
	_prev_q = _q;
	_q = new_q;
	_prev_measurement_t = _measurement_t;
	_measurement_t = new_measurement_t;
}

VectorXd model1::get_h() {

	return _h;
}

Quaterniond model1::get_q() {

	return _q;
}

Quaterniond model1::get_prev_q() {

	return _prev_q;
}

state1::state1() : state_predict() {
}

MatrixXd state1::predict_sigma(const MatrixXd& augmented_sigma, double dt) {

	MatrixXd predicted_sigma = MatrixXd(NX, _NSIGMA);

	for (int c = 0; c < _NSIGMA; ++c) {
		/*
			1. Get the current state
		*/
		const double px = augmented_sigma(0, c);
		const double py = augmented_sigma(1, c);
		const double pz = augmented_sigma(2, c);

		const double proll = augmented_sigma(3, c);
		const double ppitch = augmented_sigma(4, c);
		const double pyaw = augmented_sigma(5, c);

		const double px_noise = augmented_sigma(6, c);
		const double py_noise = augmented_sigma(7, c);
		const double pz_noise = augmented_sigma(8, c);

		const double proll_noise = augmented_sigma(9, c);
		const double ppitch_noise = augmented_sigma(10, c);
		const double pyaw_noise = augmented_sigma(11, c);

		/*
			2. Predict next state with noise
		*/
		Quaterniond old_q = AngleAxisd(fmodl(proll, 2 * M_PI), Vector3d::UnitX()) * AngleAxisd(fmodl(ppitch, 2 * M_PI), Vector3d::UnitY()) * AngleAxisd(fmodl(pyaw, 2 * M_PI), Vector3d::UnitZ());
		old_q = old_q.normalized();
		
		Eigen::MatrixXd J_e = MatrixXd(3, 4);
		
		J_e << -((model1*)_m)->get_q().x(), ((model1*)_m)->get_q().w(), -((model1*)_m)->get_q().z(),
			((model1*)_m)->get_q().y(), -((model1*)_m)->get_q().y(), ((model1*)_m)->get_q().z(), ((model1*)_m)->get_q().w(),
			-((model1*)_m)->get_q().x(), -((model1*)_m)->get_q().z(), -((model1*)_m)->get_q().y(), ((model1*)_m)->get_q().x(), ((model1*)_m)->get_q().w();
		
		Vector4d qdiff = Eigen::Vector4d(((model1*)_m)->get_q().w() - ((model1*)_m)->get_prev_q().w(),
			((model1*)_m)->get_q().x() - ((model1*)_m)->get_prev_q().x(),
			((model1*)_m)->get_q().y() - ((model1*)_m)->get_prev_q().y(),
			((model1*)_m)->get_q().z() - ((model1*)_m)->get_prev_q().z());

		Vector3d vel_rollpitchyaw = vel_rollpitchyaw = 2 * J_e * qdiff;

		double m_dt = ((model1*)_m)->get_measurement_dt();
		if (!m_dt) m_dt = 1.0;

		Vector3d w(vel_rollpitchyaw[0] / m_dt, vel_rollpitchyaw[1] / m_dt, vel_rollpitchyaw[2] / m_dt);
		Vector3d pos(px, py, pz);
		Vector3d u = pos - ((model1*)_m)->get_h();

		Vector3d vxyz = w.cross(u);
		Vector3d axyz = w.cross(w.cross(u));

		Quaterniond d_q = (old_q * Quaterniond(0, vel_rollpitchyaw[0], vel_rollpitchyaw[1], vel_rollpitchyaw[2]));
		d_q.coeffs() *= 0.5;

		Quaterniond new_q;
		new_q.coeffs() = old_q.coeffs() + d_q.coeffs() * dt;
		new_q = new_q.normalized();
		auto new_orien = new_q.toRotationMatrix().eulerAngles(0, 1, 2);
		Vector3d new_pos = (AngleAxisd(w.norm() * dt, w.normalized()).toRotationMatrix() * u) + ((model1*)_m)->get_h();

		predicted_sigma(0, c) = new_pos[0] + px_noise;
		predicted_sigma(1, c) = new_pos[1] + py_noise;
		predicted_sigma(2, c) = new_pos[2] + pz_noise;
		
		predicted_sigma(3, c) = fmodl(new_orien[0] + proll_noise, 2 * M_PI);
		predicted_sigma(4, c) = fmodl(new_orien[1] + ppitch_noise, 2 * M_PI);
		predicted_sigma(5, c) = fmodl(new_orien[2] + pyaw_noise, 2 * M_PI);
	}

	return predicted_sigma;
}

ukf1::ukf1(ALL_DATA, double t, model* m, bool debug) : ukf() {
	using namespace m1_sets;
	using namespace std;
	vector<double> noises{ VAR_PX_NOISE, VAR_PY_NOISE, VAR_PZ_NOISE, VAR_PROLL_NOISE, VAR_PPITCH_NOISE, VAR_PYAW_NOISE };
	_P = MatrixXd::Identity(NX, NX) * 1.0;

	initialize(DATA, t, NSIGMA, NAUG, W, W0_m, W0_c, noises, SCALE, VAR, m, debug);
}

state_predict& ukf1::get_stater() {
	return _stater1;
}
