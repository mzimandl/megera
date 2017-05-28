#ifndef _UNITS
	#define _UNITS

#include <cmath>

//------------------------------- SI JEDNOTKY --------------------------

const float MO = 1.98855e+30;		//kg
const float kpc = 3.086e+19;		//m
const float myr = 365.2524e+6*24*60*60;	//s
const float G = 6.672e-11;		//m^3 kg^-1 s^-2

//------------------------------ NOVE JEDNOTKY (G = 1) -----------------

const float RUNIT = 1;				//kpc - jednotka vzdálenosti
const float TUNIT = 1;				//myr - jednotka času
	const float DTMyr = 1;				//myr - integrační krok

const float VUNITms = (RUNIT*kpc)/(TUNIT*myr);	//m/s - jednotka rychlosti
const float VUNIT = VUNITms/1000;		//km/s - jednotka rychlosti

const float G2 = G/(pow(VUNITms,2)*(RUNIT*kpc));//gravitační konstanta v (kpc^3 kg^-1 myr^-2)
const double GMUNIT = 1/(G2*MO);		//hmotnost jednotka v MO, jednotky s G=1

#endif
