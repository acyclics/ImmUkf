#include "data.h"

armordata::armordata() {
}

void armordata::update_measurements(sensordata sensordata, double new_mt) {

	_prev_mt = _mt;
	_mt = new_mt;
	double m_dt = _mt - _prev_mt;
	if (!m_dt) m_dt = 1.0;

	_prev_px = _px;
	_prev_py = _py;
	_prev_pz = _pz;

	_prev_vx = _vx;
	_prev_vy = _vy;
	_prev_vz = _vz;

	_prev_ax = _ax;
	_prev_ay = _ay;
	_prev_az = _az;

	_prev_q = _q;

	_px = sensordata.px;
	_py = sensordata.py;
	_pz = sensordata.pz;

	_vx = (_px - _prev_px) / m_dt;
	_vy = (_py - _prev_py) / m_dt;
	_vz = (_pz - _prev_pz) / m_dt;

	_ax = (_vx - _prev_vx) / m_dt;
	_ay = (_vy - _prev_vy) / m_dt;
	_az = (_vz - _prev_vz) / m_dt;

	_q = Quaterniond(sensordata.qw, sensordata.qx, sensordata.qy, sensordata.qz);
}

VectorXd armordata::get_extra_measures() {
	Vector3d extra_measures = Vector3d(_vx, _vy, _vz);
	return extra_measures;
}
