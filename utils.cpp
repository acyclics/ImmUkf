#include "utils.h"
#include <cmath>

MatrixXd pdefinite_svd(MatrixXd m) {

	BDCSVD<MatrixXd> USV(m, ComputeFullU | ComputeFullV);
	double eps = 1e-11;
	VectorXd s = USV.singularValues();
	int n = USV.matrixU().rows(), p = USV.matrixV().rows();
	MatrixXd S = MatrixXd::Zero(n, p);
	for (int i = 0; i < s.size(); ++i) {
		if (s(i) <= eps)
			S(i, i) = 1e-8;
		else
			S(i, i) = s(i);
	}

	MatrixXd U = (USV.matrixU().array() < 0).select(1e-8, USV.matrixU());
	MatrixXd V = (USV.matrixV().array() < 0).select(1e-8, USV.matrixV());

	return U * S * V.transpose();
}
