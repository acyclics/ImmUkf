#include "armorfilter.h"
#include <iostream>

armorfilter::armorfilter() {
}

tda& armorfilter::get_tda() {

	return _tda;
}

void armorfilter::predict() {
	
	_tda.trackers_predict();
}

void armorfilter::predict(const double t) {

	_tda.trackers_predict(t);
}

void armorfilter::update(vector<sensordata> sensordata, vector<double> vtime) {

	MatrixXd data = MatrixXd(7, sensordata.size());

	for (int i = 0; i < sensordata.size(); ++i) {
		data(0, i) = sensordata[i].px;
		data(1, i) = sensordata[i].py;
		data(2, i) = sensordata[i].pz;
		data(3, i) = sensordata[i].qw;
		data(4, i) = sensordata[i].qx;
		data(5, i) = sensordata[i].qy;
		data(6, i) = sensordata[i].qz;
	}

	_tda.associate(data, _dist2center, vtime);
}
