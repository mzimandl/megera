#include "cObjects.h"

S_Point* cObjects::sg_object(int N) {
	if (N >= 0 && N < sg_size) {
		if (N<star.size()) return &star[N];
		else if (N<star.size()+sticky.size()) return &sticky[N-star.size()];
	}
};

S_Point* cObjects::sg_object(int G, int N) {
	if (G >= 0 && G < group.size()) if (N >= 0 && N < group[G].NStot) {
		if (N < group[G].NStars) return &star[group[G].LStars+N];
		else return &sticky[group[G].LSticky+N-group[G].NStars];
	}
};

cObjects::cObjects() {
	AddGroup();
	time = 0;
	sg_size = 0;
};

void cObjects::Reset() {
	group.clear();
	star.clear();
	sticky.clear();
	potential.clear();
	
	AddGroup();
	time = 0;
	sg_size = 0;
};

void cObjects::AddGroup() {
	S_Group temp;
	
	temp.LStars = star.size();
	temp.LSticky = sticky.size();
	temp.LPots = potential.size();
	temp.NStars = 0;
	temp.NSticky = 0;
	temp.NStot = 0;
	temp.NPots = 0;
	
	temp.relax_time = 0;
		
	group.push_back(temp);
};

void cObjects::AddStar(float gmass, float x, float y, float z, float vx, float vy, float vz) {
	S_Point temp;
	
	if (!group.empty()) {
		temp.gmass = gmass;
		temp.pos[X] = x; temp.pos[Y] = y; temp.pos[Z] = z;
		temp.vel[X] = vx; temp.vel[Y] = vy; temp.vel[Z] = vz;
			
		star.push_back(temp);
		group.back().NStars ++;
		group.back().NStot ++;
		sg_size ++;
	}
};

void cObjects::AddSticky(float gmass, float R, float x, float y, float z, float vx, float vy, float vz) {
	S_Ball temp;
	
	if (!group.empty()) {
		temp.gmass = gmass;
		temp.R = R;
		temp.pos[X] = x; temp.pos[Y] = y; temp.pos[Z] = z;
		temp.vel[X] = vx; temp.vel[Y] = vy; temp.vel[Z] = vz;
			
		sticky.push_back(temp);
		group.back().NSticky ++;
		group.back().NStot ++;
		sg_size ++;
	}
};

void cObjects::AddPotential(int type, float gmass, float chr, float x, float y, float z, float vx, float vy, float vz) {
	S_Potential temp;
	
	if (!group.empty()) {
		temp.gmass = gmass;
		temp.type = type;
		temp.chr2 = chr*chr;
		temp.pos[X] = x; temp.pos[Y] = y; temp.pos[Z] = z;
		temp.vel[X] = vx; temp.vel[Y] = vy; temp.vel[Z] = vz;
		temp.active = true;
			
		potential.push_back(temp);
		group.back().NPots ++;
	}
};

void cObjects::TurnOffPotential(int N) {
	if (N >= 0 && N < potential.size()) potential[N].active = false;
};

void cObjects::TurnOffPotential() {
	potential.back().active = false;
};

void cObjects::TurnOnPotential(int N) {
	if (N >= 0 && N < potential.size()) potential[N].active = true;
};

void cObjects::TurnOnPotential() {
	potential.back().active = true;
};

void cObjects::SetZeroPoint(float x, float y, float z) {
	int i;
	
	for (i=0;i<star.size();i++) {
		star[i].pos[X] -= x;	star[i].pos[Y] -= y;	star[i].pos[Z] -= z;
	}
	
	for (i=0;i<sticky.size();i++) {
		sticky[i].pos[X] -= x;	sticky[i].pos[Y] -= y;	sticky[i].pos[Z] -= z;
	}
	
	for (i=0;i<potential.size();i++) {
		potential[i].pos[X] -= x;	potential[i].pos[Y] -= y;	potential[i].pos[Z] -= z;
	}
};

void cObjects::SetGroupOffset(int N, float x, float y, float z, float vx, float vy, float vz) {
	int i;
	
	if (N >= 0 && N < group.size()) {
		for (i=group[N].LStars;i<group[N].LStars+group[N].NStars;i++) {
			star[i].pos[X] += x;	star[i].pos[Y] += y;	star[i].pos[Z] += z;
			star[i].vel[X] += vx;	star[i].vel[Y] += vy;	star[i].vel[Z] += vz;
		}
		
		for (i=group[N].LSticky;i<group[N].LSticky+group.back().NSticky;i++) {
			sticky[i].pos[X] += x;	sticky[i].pos[Y] += y;	sticky[i].pos[Z] += z;
			sticky[i].vel[X] += vx;	sticky[i].vel[Y] += vy;	sticky[i].vel[Z] += vz;
		}
		
		for (i=group[N].LPots;i<group[N].LPots+group[N].NPots;i++) {
			potential[i].pos[X] += x;	potential[i].pos[Y] += y;	potential[i].pos[Z] += z;
			potential[i].vel[X] += vx;	potential[i].vel[Y] += vy;	potential[i].vel[Z] += vz;
		}
	}
};

void cObjects::SetGroupOffset(float x, float y, float z, float vx, float vy, float vz) {
	int i;
	
	SetGroupOffset(group.size()-1, x, y, z, vx, vy, vz);
};

void cObjects::SetGroupRot(int N, float angle, int axis) {
	int i;
	float ta, tb;
	
	angle = 2*M_PI*angle/360;
	
	if (N >= 0 && N < group.size()) {
		switch (axis) {
			case X:
				for (i=group[N].LStars;i<group[N].LStars+group[N].NStars;i++) {
					ta = star[i].pos[Y]*cos(angle) - star[i].pos[Z]*sin(angle);
					tb = star[i].pos[Y]*sin(angle) + star[i].pos[Z]*cos(angle);
					star[i].pos[Y] = ta;	star[i].pos[Z] = tb;
					
					ta = star[i].vel[Y]*cos(angle) - star[i].vel[Z]*sin(angle);
					tb = star[i].vel[Y]*sin(angle) + star[i].vel[Z]*cos(angle);
					star[i].vel[Y] = ta;	star[i].vel[Z] = tb;
				}
		
				for (i=group[N].LSticky;i<group[N].LSticky+group.back().NSticky;i++) {
					ta = sticky[i].pos[Y]*cos(angle) - sticky[i].pos[Z]*sin(angle);
					tb = sticky[i].pos[Y]*sin(angle) + sticky[i].pos[Z]*cos(angle);
					sticky[i].pos[Y] = ta;	sticky[i].pos[Z] = tb;
					
					ta = sticky[i].vel[Y]*cos(angle) - sticky[i].vel[Z]*sin(angle);
					tb = sticky[i].vel[Y]*sin(angle) + sticky[i].vel[Z]*cos(angle);
					sticky[i].vel[Y] = ta;	sticky[i].vel[Z] = tb;
				}
		
				for (i=group[N].LPots;i<group[N].LPots+group[N].NPots;i++) {
					ta = potential[i].pos[Y]*cos(angle) - potential[i].pos[Z]*sin(angle);
					tb = potential[i].pos[Y]*sin(angle) + potential[i].pos[Z]*cos(angle);
					potential[i].pos[Y] = ta;	potential[i].pos[Z] = tb;
					
					ta = potential[i].vel[Y]*cos(angle) - potential[i].vel[Z]*sin(angle);
					tb = potential[i].vel[Y]*sin(angle) + potential[i].vel[Z]*cos(angle);
					potential[i].vel[Y] = ta;	potential[i].vel[Z] = tb;
				}
			break;
			case Y:
				for (i=group[N].LStars;i<group[N].LStars+group[N].NStars;i++) {
					ta = star[i].pos[X]*cos(angle) - star[i].pos[Z]*sin(angle);
					tb = star[i].pos[X]*sin(angle) + star[i].pos[Z]*cos(angle);
					star[i].pos[X] = ta;	star[i].pos[Z] = tb;
					
					ta = star[i].vel[X]*cos(angle) - star[i].vel[Z]*sin(angle);
					tb = star[i].vel[X]*sin(angle) + star[i].vel[Z]*cos(angle);
					star[i].vel[X] = ta;	star[i].vel[Z] = tb;
				}
		
				for (i=group[N].LSticky;i<group[N].LSticky+group.back().NSticky;i++) {
					ta = sticky[i].pos[X]*cos(angle) - sticky[i].pos[Z]*sin(angle);
					tb = sticky[i].pos[X]*sin(angle) + sticky[i].pos[Z]*cos(angle);
					sticky[i].pos[X] = ta;	sticky[i].pos[Z] = tb;
					
					ta = sticky[i].vel[X]*cos(angle) - sticky[i].vel[Z]*sin(angle);
					tb = sticky[i].vel[X]*sin(angle) + sticky[i].vel[Z]*cos(angle);
					sticky[i].vel[X] = ta;	sticky[i].vel[Z] = tb;
				}
		
				for (i=group[N].LPots;i<group[N].LPots+group[N].NPots;i++) {
					ta = potential[i].pos[X]*cos(angle) - potential[i].pos[Z]*sin(angle);
					tb = potential[i].pos[X]*sin(angle) + potential[i].pos[Z]*cos(angle);
					potential[i].pos[X] = ta;	potential[i].pos[Z] = tb;
					
					ta = potential[i].vel[X]*cos(angle) - potential[i].vel[Z]*sin(angle);
					tb = potential[i].vel[X]*sin(angle) + potential[i].vel[Z]*cos(angle);
					potential[i].vel[X] = ta;	potential[i].vel[Z] = tb;
				}
			break;
			case Z:
				for (i=group[N].LStars;i<group[N].LStars+group[N].NStars;i++) {
					ta = star[i].pos[X]*cos(angle) - star[i].pos[Y]*sin(angle);
					tb = -star[i].pos[X]*sin(angle) + star[i].pos[Y]*cos(angle);
					star[i].pos[X] = ta;	star[i].pos[Y] = tb;
					
					ta = star[i].vel[X]*cos(angle) - star[i].vel[Y]*sin(angle);
					tb = star[i].vel[X]*sin(angle) + star[i].vel[Y]*cos(angle);
					star[i].vel[X] = ta;	star[i].vel[Y] = tb;
				}
		
				for (i=group[N].LSticky;i<group[N].LSticky+group.back().NSticky;i++) {
					ta = sticky[i].pos[X]*cos(angle) - sticky[i].pos[Y]*sin(angle);
					tb = -sticky[i].pos[X]*sin(angle) + sticky[i].pos[Y]*cos(angle);
					sticky[i].pos[X] = ta;	sticky[i].pos[Y] = tb;
					
					ta = sticky[i].vel[X]*cos(angle) - sticky[i].vel[Y]*sin(angle);
					tb = sticky[i].vel[X]*sin(angle) + sticky[i].vel[Y]*cos(angle);
					sticky[i].vel[X] = ta;	sticky[i].vel[Y] = tb;
				}
		
				for (i=group[N].LPots;i<group[N].LPots+group[N].NPots;i++) {
					ta = potential[i].pos[X]*cos(angle) - potential[i].pos[Y]*sin(angle);
					tb = -potential[i].pos[X]*sin(angle) + potential[i].pos[Y]*cos(angle);
					potential[i].pos[X] = ta;	potential[i].pos[Y] = tb;
					
					ta = potential[i].vel[X]*cos(angle) - potential[i].vel[Y]*sin(angle);
					tb = potential[i].vel[X]*sin(angle) + potential[i].vel[Y]*cos(angle);
					potential[i].vel[X] = ta;	potential[i].vel[Y] = tb;
				}
			break;
		}

	}
};

void cObjects::SetGroupRot(float angle, int axis) {
	SetGroupRot(group.size()-1, angle, axis);
};

void cObjects::SetPotentialVel(int N, float vx, float vy, float vz) {
	if (N>=0 && N<potential.size()) {
		potential[N].vel[X] = vx;	potential[N].vel[Y] = vy;	potential[N].vel[Z] = vz;
	}
};