#include "ukf_plot.h"
#include <fstream>

using namespace std;

void write_residueVsTime(vector<double> time, vector<double> residue, vector<double> var, int size) {
	ofstream plot;
	plot.open("./analysis/residueVsTime");
	for (int i = 0; i < size; ++i) {
		plot << residue[i] << "," << time[i] << "," << var[i] << "\n";
	}
	plot.close();
}
