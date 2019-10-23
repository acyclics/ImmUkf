#ifndef MODEL_H
#define MODEL_H

class model {

	protected:
	double _prev_measurement_t = 0, _measurement_t = 0;

	public:
	model();
	double get_measurement_dt();
};

#endif
