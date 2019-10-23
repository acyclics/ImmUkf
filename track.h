#ifndef TRACK_H
#define TRACK_H

#include "imm.h"
#include "data.h"

class tracker {

	imm _imm;
	armordata _armordata;

	public:
	tracker();
	void initialize_imm(ALL_DATA, const VectorXd& extra_measures, double t);
	imm& get_imm();
	armordata& get_armordata();
};

#endif
