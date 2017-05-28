#ifndef _DEF_GENERATEOBJ
	#define _DEF_GENERATEOBJ

#include "math_phys.h"
#include "cObjects.h"

void GenUniformStarBox(cObjects *obj, int N, float gmass, float x1, float y1, float z1, float x2, float y2, float z2, float max_vel);
void GenUniformStickyBox(cObjects *obj, int N, float gmass, float stickyR, float x1, float y1, float z1, float x2, float y2, float z2, float max_vel);
void GenStarBall(cObjects *obj, int N, float gmass, float R, float max_vel);
void GenUniformStarBall(cObjects *obj, int N, float gmass, float R, float max_vel); // generování metodou monte-carlo
void GenHomStarBall(cObjects *obj, int N, float gmass, float R, float max_vel); // stejné jako předchozí, jiná metoda
void GenGaussStarBall(cObjects *obj, int N, float gmass, float R, float max_vel);  // koule s normálním rozdělením

void GenPlummerStar(cObjects *obj, int N, float Plum_gmass, float charR); // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
void GenPlummerSticky(cObjects *obj, int N, float Plum_gmass, float charR,  float sticky_gmass, float stickyR); // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
void GenPlummerStarTrim(cObjects *obj, int N, float Plum_gmass, float charR, float trim); // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
void GenPlummerStickyTrim(cObjects *obj, int N, float Plum_gmass, float charR,  float sticky_gmass, float stickyR, float trim); // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
void GenPlummerStarTrim2(cObjects *obj, int N, float Plum_gmass, float charR, float trim); // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
void GenPlummerStickyTrim2(cObjects *obj, int N, float Plum_gmass, float charR,  float sticky_gmass, float stickyR, float trim);

void GenPlummHommStarDisc(cObjects *obj, int N, float Plum_gmass, float charR, float R, float d);
void GenPlummHommStickyDisc(cObjects *obj, int N, float Plum_gmass, float charR, float R, float d, float sticky_gmass, float stickyR);

void GenPlummHommStarDiscHot(cObjects *obj, int N, float Plum_gmass, float charR, float R, float d);
void GenPlummHommStickyDiscHot(cObjects *obj, int N, float Plum_gmass, float charR, float R, float d, float sticky_gmass, float stickyR);

void GenPlummerStarDisc(cObjects *obj, int N, float Plum_gmass, float charR); // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
void GenPlummerStickyDisc(cObjects *obj, int N, float Plum_gmass, float charR, float stickyR); // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass

void GenStarLine(cObjects *obj, float gmass, float seg_x, float seg_y, float seg_z, int stars);

void GenDoublePlummer(cObjects *obj, float gmass1, float gmass2, float a1, float a2, int N1, int N2);
void GenDoublePlummer(cObjects *obj, float gmass1, float gmass2, float a1, float a2, int N1, int N2, float trimm1, float trimm2);

#endif