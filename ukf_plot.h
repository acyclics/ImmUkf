#ifndef UKF_PLOT_H
#define UKF_PLOT_H

#include <vector>

using namespace std;

void write_residueVsTime(vector<double> time, vector<double> residue, vector<double> var, int size);

#endif
