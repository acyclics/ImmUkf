#ifndef ARMORFILTER_H
#define ARMORFILTER_H

#include <Eigen/Dense>
#include "tda.h"
#include "global_settings.h"
#include "data.h"

using namespace std;
using namespace Eigen;

class armorfilter {

	tda _tda;
	double _dist2center = 0.5;

	public:
	armorfilter();
	tda& get_tda();
	void predict();
	void predict(const double t);
	void update(vector<sensordata> sensordata, vector<double> vtime);
};

#endif
