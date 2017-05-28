#include "AnalyzeObj.h"

void AngMomentum(cObjects &obj, float *vect) {
	vect[X]=0;
	vect[Y]=0;
	vect[Z]=0;
	
	for (int i=0;i<obj.sg_size;i++) {
		vect[X]+=obj.sg_object(i)->gmass*(obj.sg_object(i)->pos[Y]*obj.sg_object(i)->vel[Z]-obj.sg_object(i)->pos[Z]*obj.sg_object(i)->vel[Y]);
		vect[Y]+=obj.sg_object(i)->gmass*(obj.sg_object(i)->pos[Z]*obj.sg_object(i)->vel[X]-obj.sg_object(i)->pos[X]*obj.sg_object(i)->vel[Z]);
		vect[Z]+=obj.sg_object(i)->gmass*(obj.sg_object(i)->pos[X]*obj.sg_object(i)->vel[Y]-obj.sg_object(i)->pos[Y]*obj.sg_object(i)->vel[X]);
	}
}

float AngMomentumMag(cObjects *obj, int start, int end, S_Point* center, float r) {
	float vect[3] = {0,0,0};
	float temp1[3],temp2[3];
	float tot_mass = 0;
	
	for (int i=start;i<end;i++) if (Rad2(center->pos,obj->sg_object(i)->pos)<r*r) {
		for (int j=0;j<3;j++) {
			temp1[j] = obj->sg_object(i)->pos[j] - center->pos[j];
			temp2[j] = obj->sg_object(i)->vel[j] - center->vel[j];
		}
		vect[X]+=obj->sg_object(i)->gmass*(temp1[Y]*temp2[Z]-temp1[Z]*temp2[Y]);
		vect[Y]+=obj->sg_object(i)->gmass*(temp1[Z]*temp2[X]-temp1[X]*temp2[Z]);
		vect[Z]+=obj->sg_object(i)->gmass*(temp1[X]*temp2[Y]-temp1[Y]*temp2[X]);
		
		tot_mass += obj->sg_object(i)->gmass;
	}
	return sqrt(Rad2(vect))/tot_mass;
}

float AngMomentumMag(cObjects *obj, int start, int end, S_Point* center) {
	float vect[3] = {0,0,0};
	float temp1[3],temp2[3];
	float tot_mass = 0;
	
	for (int i=start;i<end;i++) {
		for (int j=0;j<3;j++) {
			temp1[j] = obj->sg_object(i)->pos[j] - center->pos[j];
			temp2[j] = obj->sg_object(i)->vel[j] - center->vel[j];
		}
		vect[X]+=obj->sg_object(i)->gmass*(temp1[Y]*temp2[Z]-temp1[Z]*temp2[Y]);
		vect[Y]+=obj->sg_object(i)->gmass*(temp1[Z]*temp2[X]-temp1[X]*temp2[Z]);
		vect[Z]+=obj->sg_object(i)->gmass*(temp1[X]*temp2[Y]-temp1[Y]*temp2[X]);
		
		tot_mass += obj->sg_object(i)->gmass;
	}
	return sqrt(Rad2(vect))/tot_mass;
}

float MeanSquareVelMag(cObjects *obj, int start, int end, S_Point* center, float r) {
	float temp = 0;
	float tot_N = 0;
	
	for (int i=start;i<end;i++) if (Rad2(center->pos,obj->sg_object(i)->pos)<r*r) {
		temp += Rad2(center->vel,obj->sg_object(i)->vel);
		tot_N++;
	}
	return temp/tot_N;
}

float MeanSquareVelMag(cObjects *obj, int start, int end, S_Point* center) {
	float temp = 0;
	float tot_N = 0;
	
	for (int i=start;i<end;i++) {
		temp += Rad2(center->vel,obj->sg_object(i)->vel);
		tot_N++;
	}
	return temp/tot_N;
}

int NParticles(cObjects *obj, int start, int end, S_Point* center, float r) {
	int tot_N = 0;
	
	for (int i=start;i<end;i++) if (Rad2(center->pos,obj->sg_object(i)->pos)<r*r) tot_N++;
	return tot_N;
}

//-------------------------

void LagrangeR(cObjects *obj, S_Point* center, std::vector<float  >& mass_frac, std::vector<float  >& mass_frac_dist) {
	std::vector<float > dists;
	
	for (int i=0; i<obj->star.size()+obj->sticky.size(); i++) dists.push_back(Rad2(center->pos,obj->sg_object(i)->pos));
	std::sort(dists.begin(), dists.end());
	
	for (int i=0; i<mass_frac.size(); i++) mass_frac_dist[i] = sqrt(dists[(int)((dists.size()-1)*mass_frac[i])]);
}

void StarLagrangeR(cObjects *obj, S_Point* center, std::vector<float  >& mass_frac, std::vector<float  >& mass_frac_dist) {
	std::vector<float > dists;
	
	for (int i=0; i<obj->star.size(); i++) dists.push_back(Rad2(center->pos,obj->star[i].pos));
	std::sort(dists.begin(), dists.end());
	
	for (int i=0; i<mass_frac.size(); i++) mass_frac_dist[i] = sqrt(dists[(int)((dists.size()-1)*mass_frac[i])]);
}

void StickyLagrangeR(cObjects *obj, S_Point* center, std::vector<float  >& mass_frac, std::vector<float  >& mass_frac_dist) {
	std::vector<float > dists;
	
	for (int i=0; i<obj->sticky.size(); i++) dists.push_back(Rad2(center->pos,obj->sticky[i].pos));
	std::sort(dists.begin(), dists.end());
	
	for (int i=0; i<mass_frac.size(); i++) mass_frac_dist[i] = sqrt(dists[(int)((dists.size()-1)*mass_frac[i])]);
}

void StarAngMomentumMag(cObjects *obj, S_Point* center, std::vector<float  >& rad, std::vector<float  >& ang_mom) {
	float vect[3];
	float temp1[3],temp2[3];
	float tot_mass;
	
	for (int j=0;j<rad.size();j++) {
		ang_mom[j] = 0;
		tot_mass = 0;
		
		vect[0] = 0;
		vect[1] = 0;
		vect[2] = 0;
		
		for (int i=0;i<obj->star.size();i++) if (Rad2(center->pos,obj->star[i].pos)<rad[j]*rad[j]) {
			for (int k=0;k<3;k++) {
				temp1[k] = obj->star[i].pos[k] - center->pos[k];
				temp2[k] = obj->star[i].vel[k] - center->vel[k];
			}
			
			vect[X]+=obj->star[i].gmass*(temp1[Y]*temp2[Z]-temp1[Z]*temp2[Y]);
			vect[Y]+=obj->star[i].gmass*(temp1[Z]*temp2[X]-temp1[X]*temp2[Z]);
			vect[Z]+=obj->star[i].gmass*(temp1[X]*temp2[Y]-temp1[Y]*temp2[X]);
		
			tot_mass += obj->star[i].gmass;
		}
		
		ang_mom[j] = sqrt(Rad2(vect))/tot_mass;
	}
}

void StickyAngMomentumMag(cObjects *obj, S_Point* center, std::vector<float  >& rad, std::vector<float  >& ang_mom) {
	float vect[3];
	float temp1[3],temp2[3];
	float tot_mass;
	
	for (int j=0;j<rad.size();j++) {
		ang_mom[j] = 0;
		tot_mass = 0;
		
		vect[0] = 0;
		vect[1] = 0;
		vect[2] = 0;
		
		for (int i=0;i<obj->sticky.size();i++) if (Rad2(center->pos,obj->sticky[i].pos)<rad[j]*rad[j]) {
			for (int k=0;k<3;k++) {
				temp1[k] = obj->sticky[i].pos[k] - center->pos[k];
				temp2[k] = obj->sticky[i].vel[k] - center->vel[k];
			}
			
			vect[X]+=obj->sticky[i].gmass*(temp1[Y]*temp2[Z]-temp1[Z]*temp2[Y]);
			vect[Y]+=obj->sticky[i].gmass*(temp1[Z]*temp2[X]-temp1[X]*temp2[Z]);
			vect[Z]+=obj->sticky[i].gmass*(temp1[X]*temp2[Y]-temp1[Y]*temp2[X]);
		
			tot_mass += obj->sticky[i].gmass;
		}
		
		ang_mom[j] = sqrt(Rad2(vect))/tot_mass;
	}
}

void StarMeanSquareVel2Mag(cObjects *obj, S_Point* center, std::vector<float  >& rad, std::vector<float  >& mean_square_vel) {
	float temp;
	float tot_N;
	
	for (int j=0;j<rad.size();j++) {
		mean_square_vel[j] = 0;
		tot_N = 0;
		
		for (int i=0;i<obj->star.size();i++) if (Rad2(center->pos,obj->star[i].pos)<rad[j]*rad[j]) {
			mean_square_vel[j] += Rad2(center->vel,obj->star[i].vel);
			tot_N++;
		}
		
		mean_square_vel[j] = mean_square_vel[j]/tot_N;
	}
	
}

void StickyMeanSquareVel2Mag(cObjects *obj, S_Point* center, std::vector<float  >& rad, std::vector<float  >& mean_square_vel) {
	float temp;
	float tot_N;
	
	for (int j=0;j<rad.size();j++) {
		mean_square_vel[j] = 0;
		tot_N = 0;
		
		for (int i=0;i<obj->sticky.size();i++) if (Rad2(center->pos,obj->sticky[i].pos)<rad[j]*rad[j]) {
			mean_square_vel[j] += Rad2(center->vel,obj->sticky[i].vel);
			tot_N++;
		}
		
		mean_square_vel[j] = mean_square_vel[j]/tot_N;
	}
	
}

//---------------------------------------

void FillDensityGrid(cObjects *obj, sGraphSetting &gs, int first, int last) {
	for (int k=first; k<last; k++) gs.data_mriz.AddData(obj->sg_object(k)->pos[X], obj->sg_object(k)->pos[Y], 1);
}

void FillPhaseGrid(cObjects *obj, sGraphSetting &gs, int first, int last, int xyz, int vxyz) {
	for (int k=first; k<last; k++) gs.data_mriz.AddData(obj->sg_object(k)->pos[xyz], obj->sg_object(k)->vel[vxyz]*VUNIT, 1);
}

void FillVelocityGrid(cObjects *obj, sGraphSetting &gs, int first, int last) {
	for (int k=first; k<last; k++) gs.data_mriz.AddData(obj->sg_object(k)->vel[X]*VUNIT, obj->sg_object(k)->vel[Y]*VUNIT, 1);
}