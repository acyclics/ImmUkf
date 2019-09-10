#ifndef M0_H
#define M0_H

#include <Eigen/Dense>
#include "state.h"
#include "ukf.h"
#include <vector>

using namespace Eigen;

class state0 : public state_predict {

	public:
	state0();
	MatrixXd predict_sigma(const MatrixXd& augmented_sigma, double dt);
};

class ukf0 : public ukf {

	state0 _stater0;

	public:
	ukf0(ALL_DATA, double t, bool debug);
	state_predict& get_stater();
};

#endif
