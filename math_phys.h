#ifndef _MATH_PHYS_DEF
	#define _MATH_PHYS_DEF

#include <cmath>
#include <cstdlib>
#include "units.h"

#define DIM 3
#define X 0
#define Y 1
#define Z 2

#define T_NEWT	0
#define T_HERN	1
#define T_PLUM	2

//----------- STRUCTURES ----------------

struct S_Point {
   float gmass;
   float pos[DIM];
   float vel[DIM];
};

struct S_Ball:S_Point {
   float R;
};

struct S_Potential:S_Point {
    int	type;
    float chr2;
    bool active;
};

struct S_StColl {
	S_Ball* st1;
	S_Ball* st2;
	float dist;
};

//--------------- MATH ------------------

float Rad2(float* x1, float* x2);
float Rad2(float* x1);

float GenUniform();
float GenUniform(float a, float b);
float GenGaussian();

//------------- PHYSICS -----------------

float AccNewt(float gmass, float x, float r2);
float AccHern(float gmass, float x, float r, float a);
float AccPlum(float gmass, float x, float r2, float a2);

float PotNewt(float gmass, float r);
float PotHern(float gmass, float r, float a);
float PotPlum(float gmass, float r, float a);
float PotPlum2(float gmass, float r2, float a2);

float VescNewt(float gmass, float r);
float VescHern(float gmass, float r, float a);
float VescPlum(float gmass, float r, float a);

float DensPlum(float gmass, float r, float a);

float DispPlum(float gmass, float r, float a);

float TidalPlum(S_Potential &pot1, S_Potential &pot2);
float InMassPlum(float r, S_Potential &pot1);
float OutMassPlum(float r, S_Potential &pot1);
float PlumMassRadius(float mass_fraction, S_Potential &pot1);

//------

float AccPotential(int type, float GM, float x, float r2, float a2);

//random radial distribution of Plummer, inverse relation of cumulative mass
float RandPlummerRad(float a);
float RandPlummerRad(float a, float trim);

//random velocity on r according to Plummer distribution function, f(x,v) ~ E^(7/2) polytrope, monte-carlo
float RandPlummerVel(float gmass, float r, float a);
float RandPlummerVel(float gmass, float r, float a, float trim);

void RandPlummerRadVel(float gmass, float a, float trim, float &r, float &v);

//handle collision, coefitients from 0 to 2 (0 - no collision, 1 - nonelastic, 2 - elastic)
void HandleCollision(S_Ball* obj1, S_Ball* obj2, float para_coef, float perp_coef);
void HandleCollision(S_StColl& balls, float para_coef, float perp_coef);

void HandleDynamicalFrictionPlummer(S_Potential &pot1, S_Potential &pot2, float *acc2);
#endif
