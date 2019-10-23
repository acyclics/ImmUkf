#ifndef IMM_H
#define IMM_H

#include "global_settings.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "ukf.h"
#include "model.h"

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
	model* _model[NM];

	void mix();
	void compute_state();
	void compute_mix_prob();
	void update_models(ALL_DATA, const VectorXd& extra_measures, double t);

	public:
	imm();
	imm(bool debug);
	~imm();
	void initialize(ALL_DATA, const VectorXd& extra_measures, double t);
	void predict(double t);
	void update(ALL_DATA, const VectorXd& extra_measures, double t);
	VectorXd peek(double dt);
	VectorXd get() const;
	VectorXd get_mu();
	VectorXd get_xi(int i);
	VectorXd get_xi_peek(int i, double dt);
	ukf& get_ukf(int i);
};

#endif
