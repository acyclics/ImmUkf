#ifndef M0_SETTINGS_H
#define M0_SETTINGS_H

#include "global_settings.h"
#include <cmath>

namespace m0_sets {
	/* Size of variables */
	const int NAUG = NX + 6;															// size of state "x" plus one for each noise variables
	const int NSIGMA = NAUG * 2 + 1;													// number of sigma points, 2N + 1

	/* State noise */
	const double STD_AX_NOISE = 0.1;
	const double STD_AY_NOISE = 0.1;
	const double STD_AZ_NOISE = 0.1;

	const double STD_PYAW_NOISE = 0.1;
	const double STD_PPITCH_NOISE = 0.1;
	const double STD_PROLL_NOISE = 0.1;

	const double VAR_AX_NOISE = STD_AX_NOISE * STD_AX_NOISE;
	const double VAR_AY_NOISE = STD_AY_NOISE * STD_AY_NOISE;
	const double VAR_AZ_NOISE = STD_AZ_NOISE * STD_AZ_NOISE;

	const double VAR_PYAW_NOISE = STD_PYAW_NOISE * STD_PYAW_NOISE;
	const double VAR_PPITCH_NOISE = STD_PPITCH_NOISE * STD_PPITCH_NOISE;
	const double VAR_PROLL_NOISE = STD_PROLL_NOISE * STD_PROLL_NOISE;

	const double alpha = 1e-3;															// tune alpha
	const double beta = 2;
	const double k = 0;																	// tune k
	const double LAMBDA = alpha * alpha * (NX + k) - NX;								// lambda
	const double SCALE = sqrt(LAMBDA + NAUG);											// used to create augmented sigma points
	const double W = 0.5 / (LAMBDA + NAUG);
	const double W0_m = LAMBDA / double(LAMBDA + NAUG);
	const double W0_c = LAMBDA / double(LAMBDA + NAUG) + (1 - alpha * alpha + beta);
}

#endif
