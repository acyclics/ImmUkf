#ifndef M2_H
#define M2_H

#include <Eigen/Dense>
#include "ukf.h"
#include "state.h"
#include <vector>

using namespace Eigen;

struct model2 {
	VectorXd h;
	Quaterniond q, prev_q;
};

class state2 : public state_predict {

	public:
	state2();
	MatrixXd predict_sigma(const MatrixXd& augmented_sigma, double dt);
};

class ukf2 : public ukf {

	state2 _stater2;

	public:
	ukf2(ALL_DATA, double t, bool debug);
	state_predict& get_stater();
};

#endif
