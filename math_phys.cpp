#include "math_phys.h"

//---------------------------- MATH ------------------------------------

//Výpočet kvadrátu vzdálenosti mezi x1 a x2(0)
float Rad2(float* x1, float* x2) { return pow(x1[X]-x2[X],2) + pow(x1[Y]-x2[Y],2) + pow(x1[Z]-x2[Z],2); }
float Rad2(float* x1) { return pow(x1[X],2) + pow(x1[Y],2) + pow(x1[Z],2); }

//Generuje float od 0 do 1, nezapomenout na srand() !!!
float GenUniform() {
	return (float)(rand())/RAND_MAX;
}

//Generuje float od a do b
float GenUniform(float a, float b) {
	return (b-a)*(float)(rand())/RAND_MAX + a;
}

//polar form Box-Muller transformation
float GenGaussian() {
	float x1,x2,w;
	
	do {
		x1 = 2*GenUniform()-1;
		x2 = 2*GenUniform()-1;
		w = x1*x1 + x2*x2;
	} while (w >= 1);
	
	w = sqrt((-2*log(w))/w);
	
	return x1*w;
	//return x2*w;
}

//-------------------------- PHYSICS -----------------------------------

//označení typů potenciálů
#define T_NEWT	0
#define T_HERN	1
#define T_PLUM	2

float AccNewt(float gmass, float x, float r2) { return -gmass*x/pow(r2,1.5); }
float AccHern(float gmass, float x, float r, float a) { return -gmass*x/pow(r+a,2); }
float AccPlum(float gmass, float x, float r2, float a2) { return -gmass*x/pow(r2+a2,1.5); }

float PotNewt(float gmass, float r) { return -gmass/r; }
float PotHern(float gmass, float r, float a) { return -gmass/(r+a); }
float PotPlum(float gmass, float r, float a) { return -gmass/sqrt(r*r+a*a); }
float PotPlum2(float gmass, float r2, float a2) { return -gmass/sqrt(r2+a2); }

float VescNewt(float gmass, float r) { return sqrt(2*gmass/r); }
float VescHern(float gmass, float r, float a) { return sqrt(2*gmass/(r+a)); }
float VescPlum(float gmass, float r, float a) { return sqrt(2*gmass/sqrt(r*r+a*a)); }

float DensPlum(float gmass, float r, float a) { return 3*gmass*a*a/(4*M_PI*pow(r*r+a*a,5/2.0)); }

float DispPlum(float gmass, float r, float a) { return sqrt(fabs(PotPlum(gmass,r,a))/6.0); }

//-----

float AccPotential(int type, float GM, float x, float r2, float a2) {
	switch (type) {
		case T_NEWT: return -GM*x/pow(r2,1.5);
		case T_HERN: return -GM*x/pow(sqrt(r2)+sqrt(a2),2);
		case T_PLUM: return -GM*x/pow(r2+a2,1.5);
	}
}

float RandPlummerRad(float a) { return a/sqrt(pow(1/GenUniform(),2./3.)-1); }
float RandPlummerVel(float gmass, float r, float a) {
	float temp1, temp2;
	do {
		temp1 = GenUniform();
		temp2 = 0.1*GenUniform();
	} while (temp2 > temp1*temp1*pow(1-temp1*temp1,3.5));
	return temp1*VescPlum(gmass,r,a);
}

void RandPlummerRadVel(float gmass, float a, float trim, float &r, float &v) {
	float temp1, temp2;
	
	do r = RandPlummerRad(a); while (r>trim);
	do {
		temp1 = GenUniform();
		temp2 = 0.1*GenUniform();
	} while (temp2 > temp1*temp1*pow(1-temp1*temp1,3.5));
	v = temp1*(VescPlum(gmass,r,a)-VescPlum(gmass,trim,a));
}

float RandPlummerRad(float a, float trim) {
	float temp;
	
	do temp = a/sqrt(pow(1/GenUniform(),2./3.)-1);
	while (temp>trim);
	
	return temp;
}
float RandPlummerVel(float gmass, float r, float a, float trim) {
	float temp1, temp2;
	float vesc = VescPlum(gmass,r,a);
	float vesc_trim = VescPlum(gmass,trim,a);
	
	do {
			temp1 = GenUniform();
			temp2 = 0.1*GenUniform();
	} while (temp2 > temp1*temp1*pow(1-temp1*temp1,3.5) || pow(temp1*vesc,2) > pow(vesc,2) - 2*pow(vesc_trim,2));
	return temp1*vesc;
}

void HandleCollision(S_Ball* obj1, S_Ball* obj2, float para_coef, float perp_coef) {
	int i;
	float temp, product, distance;
	float rel_pos[DIM], rel_vel[DIM];
	float para_vel[DIM], perp_vel[DIM];
	
	for (i=X;i<DIM;i++) rel_pos[i] = obj2->pos[i] - obj1->pos[i];
	distance = sqrt(Rad2(rel_pos));
	if (distance <= obj1->R + obj2->R && distance > 0) {
		for (i=X;i<DIM;i++) {
			rel_pos[i] /= distance;
			rel_vel[i] = obj2->vel[i] - obj1->vel[i];
		}
		
		product = 0;
		for (i=X;i<DIM;i++) product+=rel_pos[i]*rel_vel[i];
			
		for (i=X;i<DIM;i++) {
			para_vel[i] = product*rel_vel[i];
			perp_vel[i] = rel_vel[i]-para_vel[i];
		}
		
		for (i=X;i<DIM;i++) {
			temp = (para_coef*para_vel[i]+perp_coef*perp_vel[i])/(obj1->gmass+obj2->gmass);
			obj1->vel[i] += obj2->gmass*temp;
			obj2->vel[i] -= obj1->gmass*temp;
		}
	}
}

void HandleCollision(S_StColl& balls, float para_coef, float perp_coef) {
	float temp;
	float rel_pos[DIM], rel_vel[DIM];
		
	for (int i=X;i<DIM;i++) {
		rel_pos[i] = (balls.st2->pos[i] - balls.st1->pos[i])/balls.dist;
		rel_vel[i] = balls.st2->vel[i] - balls.st1->vel[i];
	}
		
	float product = 0;
	for (int i=X;i<DIM;i++) product += rel_pos[i]*rel_vel[i];
	
	float para_vel[DIM], perp_vel[DIM];		
	for (int i=X;i<DIM;i++) {
		para_vel[i] = product*rel_vel[i];	//projekce rychlosti do vzájemné spojnice
		perp_vel[i] = rel_vel[i]-para_vel[i];
		
		temp = ((1-para_coef)*para_vel[i]+(1-perp_coef)*perp_vel[i])/(balls.st1->gmass+balls.st2->gmass);
		
		balls.st1->vel[i] += balls.st2->gmass*temp;
		balls.st2->vel[i] -= balls.st1->gmass*temp;
	}
}

void HandleDynamicalFrictionPlummer(S_Potential &pot1, S_Potential &pot2, float *acc2) { //tření 2 v poli 1
	float rel_vel[3];
	for (int i=0;i<3;i++) rel_vel[i]=pot2.vel[i]-pot1.vel[i];

	float vel = sqrt(Rad2(pot1.vel,pot2.vel));
	float rad = sqrt(Rad2(pot1.pos,pot2.pos));
	
/*	float IX = 1/(sqrt(2.0)*DispPlum(pot1.gmass,sqrt(pot1.chr2),rad));
	float temp0 = 4*pow(IX,3);
	
	IX *= vel;
	float temp1 = 0;
	int factorial = 1;
	
	for (int i=1;i<100;i++) {
		factorial *= i;
		temp1 += i*pow(-1,i+1)*pow(IX,2*(i-1))/((1+2*i)*factorial);
	}

	float temp = temp0*temp1*4*sqrt(M_PI)*pot2.gmass*log(1.3*sqrt(pot1.chr2)/(0.64*sqrt(pot2.chr2)))*DensPlum(pot1.gmass,sqrt(pot1.chr2),rad); // CL je half mass radius primaru/core radius sekundaru
*/
	float IX = vel/(sqrt(2.0)*DispPlum(pot1.gmass,sqrt(pot1.chr2),rad));
	float temp = erf(IX) - 2*IX*exp(-IX*IX)/sqrt(M_PI);
	temp *= 4*M_PI*pot2.gmass*log(1.3*sqrt(pot1.chr2)/(0.64*sqrt(pot2.chr2)))*DensPlum(pot1.gmass,sqrt(pot1.chr2),rad)/pow(vel,3);
	
	for (int i=0;i<3;i++) acc2[i] += -temp*rel_vel[i];
}

float TidalPlum(S_Potential &pot1, S_Potential &pot2) { //ZHRUBA radius pot2 v poli pot1
	float rad = sqrt(Rad2(pot1.pos,pot2.pos));
	
	return sqrt((rad*rad+pot1.chr2)*pow(pot2.gmass/pot1.gmass,2/3.0)-pot2.chr2);
}

float InMassPlum(float r, S_Potential &pot1) {
	return pot1.gmass*pow(r,3)/pow(r*r+pot1.chr2,3/2.0);
}

float OutMassPlum(float r, S_Potential &pot1) {
	return pot1.gmass-InMassPlum(r,pot1);
}

float PlumMassRadius(float mass_fraction, S_Potential &pot1) {
	return sqrt(pot1.chr2/(1.0-pow(2.0/3.0,1.0/mass_fraction)));
}