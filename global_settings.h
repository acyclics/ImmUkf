#ifndef GLOBAL_SETTINGS_H
#define GLOBAL_SETTINGS_H

#define ALL_DATA const double px, const double py, const double pz, const double qw, const double qx, const double qy, const double qz

#define DATA px, py, pz, qw, qx, qy, qz

#define ALL_VAR double VAR_PX, double VAR_PY, double VAR_PZ, double VAR_PROLL, double VAR_PPITCH, double VAR_PYAW

#define VAR VAR_PX, VAR_PY, VAR_PZ, VAR_PROLL, VAR_PPITCH, VAR_PYAW

/* imm variables */
#define NM 2																// number of models for imm
#define NX 6																// size of state "x"
#define NZ 6																// size of measurement "z"

/* Measurement noise */
const double STD_PX = 0.0;
const double STD_PY = 0.0;
const double STD_PZ = 0.0;

const double STD_PROLL = 0.0;
const double STD_PPITCH = 0.0;
const double STD_PYAW = 0.0;

const double VAR_PX = STD_PX * STD_PX;
const double VAR_PY = STD_PY * STD_PY;
const double VAR_PZ = STD_PZ * STD_PZ;

const double VAR_PROLL = STD_PROLL * STD_PROLL;
const double VAR_PPITCH = STD_PPITCH * STD_PPITCH;
const double VAR_PYAW = STD_PYAW * STD_PYAW;

#endif
