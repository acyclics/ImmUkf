#ifndef M0_SETTINGS_H
#define M0_SETTINGS_H

#include "global_settings.h"
#include <cmath>

namespace m0_sets {
	/* Size of variables */
	const int NAUG = NX + 7;															// size of state "x" plus one for each noise variables
	const int NSIGMA = NAUG * 2 + 1;													// number of sigma points, 2N + 1

	/* State noise */
	const double STD_AX_NOISE = 0.9;
	const double STD_AY_NOISE = 0.9;
	const double STD_AZ_NOISE = 0.9;
	const double STD_AYAW_NOISE = 0.6;
	const double STD_APITCH_NOISE = 0.6;
	const double STD_AROLL_NOISE = 0.6;
	const double STD_DIST2CENTER_NOISE = 1;

	const double VAR_AX_NOISE = STD_AX_NOISE * STD_AX_NOISE;
	const double VAR_AY_NOISE = STD_AY_NOISE * STD_AY_NOISE;
	const double VAR_AZ_NOISE = STD_AZ_NOISE * STD_AZ_NOISE;
	const double VAR_AYAW_NOISE = STD_AYAW_NOISE * STD_AYAW_NOISE;
	const double VAR_APITCH_NOISE = STD_APITCH_NOISE * STD_APITCH_NOISE;
	const double VAR_AROLL_NOISE = STD_AROLL_NOISE * STD_AROLL_NOISE;
	const double VAR_DIST2CENTER_NOISE = STD_DIST2CENTER_NOISE * STD_DIST2CENTER_NOISE;

	/* Measurement noise */
	const double STD_PX = 0.1;
	const double STD_PY = 0.1;
	const double STD_PZ = 0.1;
	const double STD_PYAW = 0.1;
	const double STD_PPITCH = 0.1;
	const double STD_PROLL = 0.1;

	const double VAR_PX = STD_PX * STD_PX;
	const double VAR_PY = STD_PY * STD_PY;
	const double VAR_PZ = STD_PZ * STD_PZ;
	const double VAR_PYAW = STD_PYAW * STD_PYAW;
	const double VAR_PPITCH = STD_PPITCH * STD_PPITCH;
	const double VAR_PROLL = STD_PROLL * STD_PROLL;


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
