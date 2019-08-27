#ifndef STATE_H
#define STATE_H

#include "global_settings.h"
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

class state_predict {

	int _model_no;
	int _NSIGMA, _NAUG;
	vector<double> _WEIGHTS_m, _WEIGHTS_c;
	vector<double> _noises;
	int _noises_size;

	double _SCALE;

	MatrixXd _sigma;			// usually, 2N + 1 is the formula for number of sigma points
	VectorXd _x = VectorXd(NX);
	MatrixXd _P = MatrixXd(NX, NX);

	MatrixXd compute_augmented_sigma(const VectorXd& current_x, const MatrixXd& current_P);
	MatrixXd predict_sigma(const MatrixXd& augmented_sigma, double dt);
	VectorXd predict_x(const MatrixXd& predicted_sigma);
	MatrixXd predict_P(const MatrixXd& predicted_sigma, const VectorXd& predicted_x);

	public:
	state_predict();
	void initialize(int model_no, int NSIGMA, int NAUG, double W, double W0_m, double W0_c, vector<double> noises, double SCALE);
    void process(VectorXd& current_x, MatrixXd& current_P, double dt);
    MatrixXd get_sigma() const;
    VectorXd getx() const;
    MatrixXd getP() const;
};

#endif