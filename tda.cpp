#include "tda.h"
#include "hungarian.h"
#include "armorfilter.h"
#include <iostream>

tda::tda() {

}

vector<tracker>& tda::get_trackers() {

	return _trackers;
}

void tda::trackers_predict() {

	for (int i = 0; i < _trackers.size(); ++i) {
		_trackers[i].get_imm().predict(current_time());
	}
}

void tda::trackers_predict(const double t) {

	for (int i = 0; i < _trackers.size(); ++i) {
		_trackers[i].get_imm().predict(t);
	}
}

void tda::associate(const MatrixXd& sensor_data, const double dist2center, const vector<double> vtime) {
	
	MatrixXd data = MatrixXd(6, sensor_data.cols());
	for (int i = 0; i < sensor_data.cols(); ++i) {
		Quaterniond q(sensor_data(3, i), sensor_data(4, i), sensor_data(5, i), sensor_data(6, i));
		q = q.normalized();
		auto rollpitchyaw = q.toRotationMatrix().eulerAngles(0, 1, 2);

		data(0, i) = sensor_data(0, i);
		data(1, i) = sensor_data(1, i);
		data(2, i) = sensor_data(2, i);
		data(3, i) = rollpitchyaw[0];
		data(4, i) = rollpitchyaw[1];
		data(5, i) = rollpitchyaw[2];
	}

	if (_trackers.size() == 0) {
		_trackers.resize(data.cols());
		for (int i = 0; i < data.cols(); ++i) {
			VectorXd extra_measures = VectorXd(4);
			extra_measures << dist2center, _trackers[i].get_armordata().get_extra_measures();
			_trackers[i].initialize_imm(sensor_data(0, i), sensor_data(1, i), sensor_data(2, i), sensor_data(3, i),
				sensor_data(4, i), sensor_data(5, i), sensor_data(6, i), extra_measures, vtime[i]);
		}
	}

	MatrixXd cost(_trackers.size(), data.cols());
	for (int i = 0; i < _trackers.size(); ++i) {
		for (int j = 0; j < data.cols(); ++j) {
			VectorXd diff = _trackers[i].get_imm().get() - data.col(j);
			double l2_cost = diff.norm();
			cost(i, j) = l2_cost;
		}
	}

	vector<int> assignments = hungarian_solve(cost);
	vector<int> unassigned_tracks;

	for (int i = 0; i < assignments.size(); ++i) {
		if (assignments[i] != -1) {
			if (cost(i, assignments[i]) > GATE_DIST) {
				// delete track
				assignments[i] = -1;
				unassigned_tracks.push_back(i);
			}
		}
	}

	// check for unassigned data
	vector<int> assigned_data(data.cols(), 0);
	for (int i = 0; i < assignments.size(); ++i) {
		assigned_data[assignments[i]] |= 1;
	}

	for (int i = 0; i < assigned_data.size(); ++i) {
		if (!assigned_data[i]) {
			_trackers.push_back(tracker());
		}
	}

	// update old tracks and initialize new ones
	for (int i = 0; i < assignments.size(); ++i) {
		sensordata sd = { sensor_data(0, i), sensor_data(1, i), sensor_data(2, i), sensor_data(3, i),
						  sensor_data(4, i), sensor_data(5, i), sensor_data(6, i) };
		_trackers[i].get_armordata().update_measurements(sd, vtime[assignments[i]]);
		
		VectorXd extra_measures = VectorXd(4);
		extra_measures << dist2center, _trackers[i].get_armordata().get_extra_measures();
		
		if (assignments[i] != -1) {
			_trackers[i].get_imm().update(sensor_data(0, assignments[i]), sensor_data(1, assignments[i]),
										  sensor_data(2, assignments[i]), sensor_data(3, assignments[i]),
										  sensor_data(4, assignments[i]), sensor_data(5, assignments[i]),
										  sensor_data(6, assignments[i]),
										  extra_measures,
										  vtime[assignments[i]]);
		}
		else {
			_trackers[i].get_imm().initialize(sensor_data(0, assignments[i]), sensor_data(1, assignments[i]),
											  sensor_data(2, assignments[i]), sensor_data(3, assignments[i]),
										      sensor_data(4, assignments[i]), sensor_data(5, assignments[i]),
											  sensor_data(6, assignments[i]),
											  extra_measures,
											  vtime[assignments[i]]);
		}
	}
}

double tda::current_time() {

	const std::chrono::duration<double> duration = std::chrono::system_clock::now() - tstart;
	return duration.count();
}
