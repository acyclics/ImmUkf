#ifndef M3_H
#define M3_H

#include <Eigen/Dense>
#include "ukf.h"
#include "state.h"
#include "model.h"
#include <vector>

using namespace Eigen;

class model3 : public model {

	Vector3d _current_vxyz = Vector3d::Zero();
	Vector3d _measure_rollpitchyaw = Vector3d::Zero();

	public:
	model3();
	void update(double vx, double vy, double vz, double qw, double qx, double qy, double qz, double new_measurement_t);
	Vector3d get_current_vxyz();
	Vector3d get_measure_rollpitchyaw();
};

class state3 : public state_predict {

public:
	state3();
	MatrixXd predict_sigma(const MatrixXd& augmented_sigma, double dt);
};

class ukf3 : public ukf {

	state3 _stater3;

	public:
	ukf3(ALL_DATA, double t, model* m, bool debug);
	state_predict& get_stater();
};

#endif
