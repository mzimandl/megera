#include "units.h"

//hmotnosti - M0/GMUNIT
#define GMP	3.2e11/GMUNIT		//hmotnost primáru v jednotkách integrátoru
#define GMStoP	0.01			//poměr hmotnosti sekundáru a primáru
#define GMS	GMP*GMStoP		//hmotnost sekundáru

//poloměry - kpc
#define PPlumSc	5e0			//Plummerův poloměr primáru
#define SPlumSc	0.5			//Plummerův poloměr sekundáru

//polohy - kpc
//#define Dini	91.2			//počáteční separace
#define Dini	100			//počáteční separace
#define DTiniP	Dini*GMStoP/(1+GMStoP)	//poloha primáru od těžiště
#define DTiniS	Dini*1/(1+GMStoP)	//poloha sekundáru od těžiště
