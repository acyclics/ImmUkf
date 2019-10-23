#include "utils.h"
#include <cmath>

MatrixXd pdefinite_svd(MatrixXd m) {

	JacobiSVD<MatrixXd> USV(m, ComputeFullU | ComputeFullV);
	double eps = 1e-11;
	VectorXd s = USV.singularValues();
	int n = USV.matrixU().rows(), p = USV.matrixV().rows();
	MatrixXd S = MatrixXd::Zero(n, p);
	for (int i = 0; i < s.size(); ++i) {
		if (s(i) < 0)
			S(i, i) = 1e-10;
		else
			S(i, i) = s(i);
	}

	//MatrixXd U = (USV.matrixU().array() < 0).select(1e-10, USV.matrixU());
	//MatrixXd V = (USV.matrixV().array() < 0).select(1e-10, USV.matrixV());

	return USV.matrixU()* S* USV.matrixV().transpose();
}

double normalize_angle(double angle) {
	angle = fmod(angle, 2 * M_PI);
	return (angle > M_PI) ? (angle - 2 * M_PI) : (angle);
}
