#ifndef TDA_H
#define TDA_H

#include <vector>
#include <Eigen/Dense>
#include "track.h"
#include <chrono>

const double GATE_DIST = 100.0;

using namespace std;
using namespace Eigen;

class tda {

	vector<tracker> _trackers;

	public:
	chrono::time_point< chrono::system_clock> tstart = std::chrono::system_clock::now();
	tda();
	vector<tracker>& get_trackers();
	void trackers_predict();
	void trackers_predict(const double t);
	void associate(const MatrixXd& sensor_data, const double dist2center, const vector<double> vtime);
	double current_time();
};

#endif
