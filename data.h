#ifndef DATA_H
#define DATA_H

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class sensordata {
public:
	double px = 0, py = 0, pz = 0, qw = 1, qx = 0, qy = 0, qz = 0;
};

class armordata {

	Quaterniond _q = Quaterniond(1, 0, 0, 0), _prev_q = Quaterniond(1, 0, 0, 0);
	double _px = 0, _py = 0, _pz = 0, _prev_px = 0, _prev_py = 0, _prev_pz = 0;
	double _vx = 0, _vy = 0, _vz = 0, _prev_vx = 0, _prev_vy = 0, _prev_vz = 0;
	double _ax = 0, _ay = 0, _az = 0, _prev_ax = 0, _prev_ay = 0, _prev_az = 0;
	double _mt = 0, _prev_mt = 0;

	public:
	armordata();
	void update_measurements(sensordata sensordata, double new_mt);
	VectorXd get_extra_measures();
};

#endif
