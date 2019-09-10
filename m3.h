#ifndef M3_H
#define M3_H

#include <Eigen/Dense>
#include "ukf.h"
#include "state.h"
#include <vector>

using namespace Eigen;

class state3 : public state_predict {

public:
	state3();
	MatrixXd predict_sigma(const MatrixXd& augmented_sigma, double dt);
};

class ukf3 : public ukf {

	state3 _stater3;

	public:
	ukf3(ALL_DATA, double t, bool debug);
	state_predict& get_stater();
};

#endif
