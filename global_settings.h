#ifndef GLOBAL_SETTINGS_H
#define GLOBAL_SETTINGS_H

#define ALL_DATA const double px, const double py, const double pz, const double yaw, const double pitch, const double roll, const double xvel, const double yvel, const double zvel, const double yawvel, const double pitchvel, const double rollvel, const double xacc, const double yacc, const double zacc, const double yawacc, const double pitchacc, const double rollacc, const double dist2center

#define DATA px, py, pz, yaw, pitch, roll, xvel, yvel, zvel, yawvel, pitchvel, rollvel, xacc, yacc, zacc, yawacc, pitchacc, rollacc, dist2center

#define ALL_VAR double VAR_PX, double VAR_PY, double VAR_PZ, double VAR_PYAW, double VAR_PPITCH, double VAR_PROLL, double VAR_XVEL, double VAR_YVEL, double VAR_ZVEL, double VAR_YAWVEL, double VAR_PITCHVEL, double VAR_ROLLVEL, double VAR_XACC, double VAR_YACC, double VAR_ZACC, double VAR_YAWACC, double VAR_PITCHACC, double VAR_ROLLACC, double VAR_DIST2CENTER

#define VAR VAR_PX, VAR_PY, VAR_PZ, VAR_PYAW, VAR_PPITCH, VAR_PROLL, VAR_XVEL, VAR_YVEL, VAR_ZVEL, VAR_YAWVEL, VAR_PITCHVEL, VAR_ROLLVEL, VAR_XACC, VAR_YACC, VAR_ZACC, VAR_YAWACC, VAR_PITCHACC, VAR_ROLLACC, VAR_DIST2CENTER

/* imm variables */
#define NM 4																// number of models for imm
#define NX 19																// size of state "x"
#define NZ 19																// size of measurement "z"

/* Measurement noise */
const double STD_PX = 0.0;
const double STD_PY = 0.0;
const double STD_PZ = 0.0;

const double STD_PYAW = 0.0;
const double STD_PPITCH = 0.0;
const double STD_PROLL = 0.0;

const double STD_XVEL = 0.01;
const double STD_YVEL = 0.01;
const double STD_ZVEL = 0.01;

const double STD_YAWVEL = 0.001;
const double STD_PITCHVEL = 0.001;
const double STD_ROLLVEL = 0.001;

const double STD_XACC = 0.0001;
const double STD_YACC = 0.0001;
const double STD_ZACC = 0.0001;

const double STD_YAWACC = 0.0001;
const double STD_PITCHACC = 0.0001;
const double STD_ROLLACC = 0.0001;

const double STD_DIST2CENTER = 0.1;

const double VAR_PX = STD_PX * STD_PX;
const double VAR_PY = STD_PY * STD_PY;
const double VAR_PZ = STD_PZ * STD_PZ;

const double VAR_PYAW = STD_PYAW * STD_PYAW;
const double VAR_PPITCH = STD_PPITCH * STD_PPITCH;
const double VAR_PROLL = STD_PROLL * STD_PROLL;

const double VAR_XVEL = STD_XVEL * STD_XVEL;
const double VAR_YVEL = STD_YVEL * STD_YVEL;
const double VAR_ZVEL = STD_ZVEL * STD_ZVEL;

const double VAR_YAWVEL = STD_YAWVEL * STD_YAWVEL;
const double VAR_PITCHVEL = STD_PITCHVEL * STD_PITCHVEL;
const double VAR_ROLLVEL = STD_ROLLVEL * STD_ROLLVEL;

const double VAR_XACC = STD_XACC * STD_XACC;
const double VAR_YACC = STD_YACC * STD_YACC;
const double VAR_ZACC = STD_ZACC * STD_ZACC;

const double VAR_YAWACC = STD_YAWACC * STD_YAWACC;
const double VAR_PITCHACC = STD_PITCHACC * STD_PITCHACC;
const double VAR_ROLLACC = STD_ROLLACC * STD_ROLLACC;

const double VAR_DIST2CENTER = STD_DIST2CENTER * STD_DIST2CENTER;

/* imm transition probabilities */

#endif
