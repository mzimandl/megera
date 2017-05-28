#ifndef _DEF_ANALYZE
	#define _DEF_ANALYZE

#include "cObjects.h"
#include "sSettings.h"
#include "units.h"

void AngMomentum(cObjects &obj, float *vect);
float AngMomentumMag(cObjects *obj, int start, int end, S_Point* center, float r);
float AngMomentumMag(cObjects *obj, int start, int end, S_Point* center);
float MeanSquareVelMag(cObjects *obj, int start, int end, S_Point* center, float r);
float MeanSquareVelMag(cObjects *obj, int start, int end, S_Point* center);
int NParticles(cObjects *obj, int start, int end, S_Point* center, float r);

#include <vector>
#include <algorithm>
void LagrangeR(cObjects *obj, S_Point* center, std::vector<float  >& mass_frac, std::vector<float  >& mass_frac_dist);
void StarLagrangeR(cObjects *obj, S_Point* center, std::vector<float  >& mass_frac, std::vector<float  >& mass_frac_dist);
void StickyLagrangeR(cObjects *obj, S_Point* center, std::vector<float  >& mass_frac, std::vector<float  >& mass_frac_dist);

void StarAngMomentumMag(cObjects *obj, S_Point* center, std::vector<float  >& mass_rad, std::vector<float  >& ang_mom);
void StarMeanSquareVel2Mag(cObjects *obj, S_Point* center, std::vector<float  >& mass_rad, std::vector<float  >& mean_square_vel);
void StickyAngMomentumMag(cObjects *obj, S_Point* center, std::vector<float  >& mass_rad, std::vector<float  >& ang_mom);
void StickyMeanSquareVel2Mag(cObjects *obj, S_Point* center, std::vector<float  >& mass_rad, std::vector<float  >& mean_square_vel);

// plnění hustotní sítě
void FillDensityGrid(cObjects *obj, sGraphSetting &gs, int first, int last);
void FillPhaseGrid(cObjects *obj, sGraphSetting &gs, int first, int last, int xyz, int vxyz);
void FillVelocityGrid(cObjects *obj, sGraphSetting &gs, int first, int last);

#endif