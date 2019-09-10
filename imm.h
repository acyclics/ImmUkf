#ifndef IMM_H
#define IMM_H

#include "global_settings.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "ukf.h"

using namespace Eigen;

class imm {
	
	bool _initialized;
	bool _debug;
	VectorXd _x;
	MatrixXd _p;
	MatrixXd _T;
	VectorXd _mu;
	VectorXd _llh;
	MatrixXd _omega;
	MatrixXd _X;
	VectorXd _normalizer;
	std::vector<Matrix<double, NX, NX>, Eigen::aligned_allocator<Matrix<double, NX, NX> > > _P;

	ukf* _ukf[NM];

	void initialize(ALL_DATA, double t);
	void mix(MatrixXd models_x);
	void compute_state();
	void compute_mix_prob();
	void compute(ALL_DATA, double t);

	public:
	imm(bool debug);
	~imm();
	void process(ALL_DATA, double t);	// process is predict + update
	void predict(double t);
	void update(ALL_DATA, double t);
	VectorXd peek(double dt);
	VectorXd get() const;
	VectorXd get_mu();
	VectorXd get_xi(int i);
	VectorXd get_xi_peek(int i, double dt);
};

#endif
