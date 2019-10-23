#include "track.h"

tracker::tracker() {
}

void tracker::initialize_imm(ALL_DATA, const VectorXd& extra_measures, double t) {

	_imm.initialize(DATA, extra_measures, t);
}

imm& tracker::get_imm() {

	return _imm;
}

armordata& tracker::get_armordata() {

	return _armordata;
}
