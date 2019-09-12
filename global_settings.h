#ifndef GLOBAL_SETTINGS_H
#define GLOBAL_SETTINGS_H

#define ALL_DATA const double px, const double py, const double pz, const double yaw, const double pitch, const double roll, const double dist2center

#define DATA px, py, pz, yaw, pitch, roll, dist2center

#define ALL_VAR double VAR_PX, double VAR_PY, double VAR_PZ, double VAR_PYAW, double VAR_PPITCH, double VAR_PROLL, double VAR_DIST2CENTER

#define VAR VAR_PX, VAR_PY, VAR_PZ, VAR_PYAW, VAR_PPITCH, VAR_PROLL, VAR_DIST2CENTER

/* imm variables */
#define NM 2																// number of models for imm
#define NX 19																// size of state "x"
#define NZ 7																// size of measurement "z"

/* Measurement noise */
const double STD_PX = 0.0;
const double STD_PY = 0.0;
const double STD_PZ = 0.0;

const double STD_PYAW = 0.0001;
const double STD_PPITCH = 0.0001;
const double STD_PROLL = 0.0001;

const double STD_DIST2CENTER = 0.1;

const double VAR_PX = STD_PX * STD_PX;
const double VAR_PY = STD_PY * STD_PY;
const double VAR_PZ = STD_PZ * STD_PZ;

const double VAR_PYAW = STD_PYAW * STD_PYAW;
const double VAR_PPITCH = STD_PPITCH * STD_PPITCH;
const double VAR_PROLL = STD_PROLL * STD_PROLL;

const double VAR_DIST2CENTER = STD_DIST2CENTER * STD_DIST2CENTER;

#endif
