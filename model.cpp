#include "model.h"

model::model() {
}

double model::get_measurement_dt() {
	return _measurement_t - _prev_measurement_t;
}
