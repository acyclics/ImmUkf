#ifndef IMM_H
#define IMM_H

#include "global_settings.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "ukf.h"

using namespace Eigen;

class imm {
	
	bool _initialized;
	VectorXd _x;
	MatrixXd _T = MatrixXd(NM, NM);
	MatrixXd _u = MatrixXd(NM, NM);
	MatrixXd _X = VectorXd(NX, NM);
	std::vector<Matrix<double, NX, NX>, Eigen::aligned_allocator<Matrix<double, NX, NX> > > _P;

	ukf _ukf[NM];

	void initialize(ALL_DATA, double t);
	void mix(MatrixXd models_x);
	void filtering(ALL_DATA, double t);
	void compute(ALL_DATA, double t);

	public:
	imm();
	void process(ALL_DATA, double t);
};

#endif
