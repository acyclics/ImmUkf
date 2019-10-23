#ifndef M1_H
#define M1_H

#include <Eigen/Dense>
#include "ukf.h"
#include "state.h"
#include "model.h"
#include <vector>

using namespace Eigen;

class model1 : public model {

	Vector3d _h = Vector3d::Zero();
	Quaterniond _q = Quaterniond(1, 0, 0, 0), _prev_q = Quaterniond(1, 0, 0, 0);

	public:
	model1();
	void update(VectorXd& new_h, Quaterniond& new_q, double new_measurement_t);
	VectorXd get_h();
	Quaterniond get_q();
	Quaterniond get_prev_q();
};

class state1 : public state_predict {

	public:
	state1();
	MatrixXd predict_sigma(const MatrixXd& augmented_sigma, double dt);
};

class ukf1 : public ukf {

	state1 _stater1;

	public:
	ukf1(ALL_DATA, double t, model* m, bool debug);
	state_predict& get_stater();
};

#endif
