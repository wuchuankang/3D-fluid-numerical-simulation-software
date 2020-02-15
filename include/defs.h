//user-defined library-------------------

#ifndef DEFS_H
#define DEFS_H

#include <cmath>
#include <cfloat>

enum class TurboTypes{Axial, Centrifugal, General};

enum class EquationTypes{Euler, NavierStokes};

enum class InitialTypes{Uniform, Interpolate};

enum class TimestepTypes{Local, Global};

enum class TimeDiscrTypes{Runge_Kutta3, Runge_Kutta4, Runge_Kutta5, Runge_Kutta5_Hybrid, LU_SGS};

enum class SpaceDiscrTypes{JST, L_F, StegerWarming, Roe, Van_Leer, AUSM_Plus, AUSM_PW, HLLC};

enum class ReconstructTypes{None, NND, MUSCL2u, MUSCL2c, MUSCL2o, MUSCL3};

enum class ReconstructVars{ConsVars, PrimVars, CharactVars};

enum class TurbulenceTypes{None, BL, SA, SST};

// floating point type (SGLPREC=double precision, otherwise single) ***********
#ifndef SINGLE
  typedef float  REAL;        /**< floating-point value */
  #define EPSGLO FLT_EPSILON  /**< tolerance for a small number */
  #define ABS    fabsf
  #define SQRT   sqrtf
  #define SIN    sinf
  #define COS    cosf
  #define EXP    expf
  #define POW    powf
  #define TAN    tanf
  #define ATAN   atanf
  #define ATAN2  atan2f
  #define LOG10  log10f
#else
  typedef double REAL;        /**< floating-point value */
  #define EPSGLO DBL_EPSILON  /**< tolerance for a small number */
  #define ABS    fabs
  #define SQRT   sqrt
  #define SIN    sin
  #define COS    cos
  #define EXP    exp
  #define POW    pow
  #define TAN    tan
  #define ATAN   atan
  #define ATAN2  atan2
  #define LOG10  log10
#endif

// general structures of geometry and flow variables **************************
//coordinates
struct TNode
{
	REAL x, y, z;
	TNode();		// constructor function
};

struct TFace
{
	REAL	r, th, s;
	TNode	coord,
			S;
	REAL	gradx_kesi, grady_kesi, gradz_kesi,
			gradx_yita, grady_yita,	gradz_yita,
			gradx_zita, grady_zita,	gradz_zita;

	TFace();
};

struct TCell
{
	REAL    dis_Wall;
	REAL	r, th, Vol;	//wall distance
	TNode	coord;
	REAL	gradx_kesi,	grady_kesi, gradz_kesi,//The gradient of kesi
			gradx_yita,	grady_yita, gradz_yita,//The gradient of yita
			gradx_zita,	grady_zita, gradz_zita;//The gradient of zita	

	TCell();
};
//conservation variables
struct TConsVars
{
	REAL dens, xmom, ymom, zmom, ener;

	TConsVars();
};

// primitive variables---which can be obsevered directly by experment.
// so,on include the internal enery
struct TPrimVars
{
	REAL dens, uvel, vvel, wvel, press, temp;

	TPrimVars();
};

// Dependent variables (inviscid flow).
struct TVisVars
{
	REAL cond, mue;

	TVisVars();
};

//turbulence eddy viscostity
struct TEddyVars
{
	REAL ev;

	TEddyVars();
};

// override operators
// about TConsVars
TConsVars operator+ (const TConsVars &A, const TConsVars &B);

TConsVars operator+ (const TConsVars &A, REAL b);

TConsVars operator+ (REAL b, const TConsVars &A);

TConsVars operator- (const TConsVars &A, const TConsVars &B);

TConsVars operator- (const TConsVars &A, REAL b);

TConsVars operator- (REAL b, const TConsVars &A);

TConsVars operator* (const TConsVars &A, const TConsVars &B);

TConsVars operator* (const TConsVars &A, REAL b);

TConsVars operator* (REAL b, const TConsVars &A);	

TConsVars operator/ (const TConsVars &A, const TConsVars &B);

TConsVars operator/ (const TConsVars &A, REAL b );

TConsVars operator/ (REAL b, const TConsVars &A );

TConsVars POW(const TConsVars &A, REAL b );

// about TPrimVars
TPrimVars operator+ (const TPrimVars &A, const TPrimVars &B);

TPrimVars operator+ (const TPrimVars &A, REAL b);

TPrimVars operator+ (REAL b, const TPrimVars &A);

TPrimVars operator- (const TPrimVars &A, const TPrimVars &B);

TPrimVars operator- (REAL b, const TPrimVars &A);

TPrimVars operator* (const TPrimVars &A, const TPrimVars &B);

TPrimVars operator* (REAL b, const TPrimVars &A);

TPrimVars operator/ (const TPrimVars &A, REAL b);

TPrimVars operator/ (REAL b, const TPrimVars &A);

TPrimVars operator/ (const TPrimVars &A, const TPrimVars &B);

TPrimVars POW(const TPrimVars &A, REAL b );

// about TNode
REAL operator* (const TNode &A, const TNode &B);

TNode operator- (const TNode &A, const TNode &B);

TNode operator+ (const TNode &A, const TNode &B);

TNode operator* (const TNode &A, REAL b);

TNode operator* (REAL b, const TNode &A);

TNode operator/ (const TNode &A, REAL b);

//matrix product
TConsVars operator* (REAL A[][6], TConsVars &b);

// max among four values
REAL MAX4(REAL a1, REAL a2, REAL a3, REAL a4);
// max among three values
REAL MAX3(REAL a1, REAL a2, REAL a3);
//
int  MAX3(int a1, int a2, int a3);
//
int MIN3(int a1, int a2, int a3);
// min among three values
REAL MIN3(REAL a1, REAL a2, REAL a3);

//void Thomas(int TMAX, REAL *a, REAL *b, REAL *c, REAL *d, REAL *x);

//number of threads
extern int NT;

//general constants ****************************
#define PI  3.1415926535897932f
#define RAD 1.7453292519943296e-2f

//function macros*******************************
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define SIGN(a,b) (((b) > (0)) ? (a) : (-a))

#endif