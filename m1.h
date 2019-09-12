#ifndef M1_H
#define M1_H

#include <Eigen/Dense>
#include "ukf.h"
#include "state.h"
#include <vector>

using namespace Eigen;

struct model1 {
	VectorXd h;
	Quaterniond q, prev_q;
};

class state1 : public state_predict {

	public:
	state1();
	MatrixXd predict_sigma(const MatrixXd& augmented_sigma, double dt);
};

class ukf1 : public ukf {

	state1 _stater1;

	public:
	ukf1(ALL_DATA, double t, bool debug);
	state_predict& get_stater();
};

#endif
