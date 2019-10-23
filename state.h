#ifndef STATE_H
#define STATE_H

#include "global_settings.h"
#include <Eigen/Dense>
#include <vector>
#include "model.h"

using namespace Eigen;
using namespace std;

class state_predict {

	protected:
	int _NSIGMA, _NAUG;
	vector<double> _WEIGHTS_m, _WEIGHTS_c;
	vector<double> _noises;
	int _noises_size;

	double _SCALE;

	MatrixXd _sigma;			// usually, 2N + 1 is the formula for number of sigma points
	VectorXd _x;
	MatrixXd _P;

	model* _m;

	MatrixXd compute_augmented_sigma(const VectorXd& current_x, const MatrixXd& current_P);
	virtual MatrixXd predict_sigma(const MatrixXd& augmented_sigma, double dt);
	VectorXd predict_x(const MatrixXd& predicted_sigma);
	MatrixXd predict_P(const MatrixXd& predicted_sigma, const VectorXd& predicted_x);

	public:
	state_predict();
	void initialize(int NSIGMA, int NAUG, double W, double W0_m, double W0_c, vector<double> noises, double SCALE, model* m);
	void process(const VectorXd& current_x, const MatrixXd& current_P, double dt);
	VectorXd peek(const VectorXd& current_x, const MatrixXd& current_P, double dt);
    MatrixXd get_sigma() const;
    VectorXd getx() const;
    MatrixXd getP() const;
};

#endif
