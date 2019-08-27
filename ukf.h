#ifndef UKF_H
#define UKF_H

#include "global_settings.h"
#include <Eigen/Dense>
#include "state.h"
#include "measurement.h"
#include "merge.h"

class ukf {

	int _model_no;
	bool _initialized;
	double _t;
	double _nis;
	VectorXd _x = VectorXd(NX);
	MatrixXd _P = MatrixXd(NX, NX);
	state_predict _stater;
	measurement_predict _measurer;
	merger _merger;

	void update(ALL_DATA, double t, VectorXd imm_x, MatrixXd imm_P);

	public:
	ukf();
	void initialize(ALL_DATA, double t, int model_no, int NSIGMA, int NAUG, double W, double W0_m, double W0_c, vector<double> noises, double SCALE,
					double VAR_PX, double VAR_PY, double VAR_PZ, double VAR_PYAW, double VAR_PPITCH, double VAR_PROLL);
	void process(ALL_DATA, double t, VectorXd imm_x, MatrixXd imm_P);
	VectorXd get() const;
	double get_nis() const;
};

#endif
