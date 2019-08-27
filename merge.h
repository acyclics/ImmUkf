#ifndef MERGE_H
#define MERGE_H

#include "global_settings.h"
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

class merger {
	
	int _NSIGMA;
	vector<double> _WEIGHTS_m, _WEIGHTS_c;

	MatrixXd _x;
	MatrixXd _P;
	double _nis;

	MatrixXd compute_Tc(const VectorXd& predicted_x, const VectorXd& predicted_z, const MatrixXd& sigma_x, const MatrixXd& sigma_z);
	void update(const VectorXd& z, const MatrixXd& S, const MatrixXd& Tc, const VectorXd& predicted_z, const VectorXd& predicted_x, const MatrixXd& predicted_P);

	public:
	merger();
	void initialize(int NSIGMA, double W, double W0_m, double W0_c);
	void process(const VectorXd& predicted_x, const VectorXd& predicted_z, const VectorXd& z, const MatrixXd& S, const MatrixXd& predicted_P,
				 const MatrixXd& sigma_x, const MatrixXd& sigma_z);
	VectorXd getx() const;
	MatrixXd getP() const;
	double get_nis() const;
};

#endif
