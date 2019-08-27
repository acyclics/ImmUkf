#include "utils.h"
#include <cmath>

double normalize(const double a) {
	return (fabs(a) > M_PI) ? remainder(a, 2. * M_PI) : a;
}
