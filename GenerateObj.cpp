#include "GenerateObj.h"

void GenUniformStarBox(cObjects *obj, int N, float gmass, float x1, float y1, float z1, float x2, float y2, float z2, float max_vel) {
	int i;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		temp.pos[X] = (x2-x1)*GenUniform() + x1;
		temp.pos[Y] = (y2-y1)*GenUniform() + y1;
		temp.pos[Z] = (z2-z1)*GenUniform() + z1;
		temp.vel[X] = max_vel*(2*GenUniform() - 1);
		temp.vel[Y] = max_vel*(2*GenUniform() - 1);
		temp.vel[Z] = max_vel*(2*GenUniform() - 1);
		
		obj->AddStar(gmass,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenUniformStickyBox(cObjects *obj, int N, float gmass, float stickyR, float x1, float y1, float z1, float x2, float y2, float z2, float max_vel) {
	int i;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		temp.pos[X] = (x2-x1)*GenUniform() + x1;
		temp.pos[Y] = (y2-y1)*GenUniform() + y1;
		temp.pos[Z] = (z2-z1)*GenUniform() + z1;
		temp.vel[X] = max_vel*(2*GenUniform() - 1);
		temp.vel[Y] = max_vel*(2*GenUniform() - 1);
		temp.vel[Z] = max_vel*(2*GenUniform() - 1);
		
		obj->AddSticky(gmass,stickyR,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenStarBall(cObjects *obj, int N, float gmass, float R, float max_vel) {
	int i;
	float pom,alfa,beta;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		pom = R*GenUniform();
		alfa = 2*M_PI*GenUniform();
		beta = M_PI*GenUniform();
		
		temp.pos[X] = pom*cos(alfa)*sin(beta);
		temp.pos[Y] = pom*sin(alfa)*sin(beta);
		temp.pos[Z] = pom*cos(beta);
		
		pom = max_vel*GenUniform();
		alfa = 2*M_PI*GenUniform();
		beta = M_PI*GenUniform();
		
		temp.vel[X] = pom*cos(alfa)*sin(beta);
		temp.vel[Y] = pom*sin(alfa)*sin(beta);
		temp.vel[Z] = pom*cos(beta);
		
		obj->AddStar(gmass,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenUniformStarBall(cObjects *obj, int N, float gmass, float R, float max_vel) {  // generování metodou monte-carlo
	int i;
	float pom,alfa,beta;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		do {
			temp.pos[X] = R*(2*GenUniform()-1);
			temp.pos[Y] = R*(2*GenUniform()-1);
			temp.pos[Z] = R*(2*GenUniform()-1);
		} while (sqrt(Rad2(temp.pos))>R);
		
		do {
			temp.vel[X] = max_vel*(2*GenUniform()-1);
			temp.vel[Y] = max_vel*(2*GenUniform()-1);
			temp.vel[Z] = max_vel*(2*GenUniform()-1);
		} while (sqrt(Rad2(temp.vel))>max_vel);
		
		obj->AddStar(gmass,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenHomStarBall(cObjects *obj, int N, float gmass, float R, float max_vel) {  // stejné jako předchozí, jiná metoda
	int i;
	float pom,alfa,sintheta,costheta;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		pom = R*pow(GenUniform(),1./3.);
		alfa = 2*M_PI*GenUniform();
		
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.pos[X] = pom*cos(alfa)*sintheta;
		temp.pos[Y] = pom*sin(alfa)*sintheta;
		temp.pos[Z] = pom*costheta;
		
		pom = max_vel*pow(GenUniform(),1./3.);
		alfa = 2*M_PI*GenUniform();
		
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.vel[X] = pom*cos(alfa)*sintheta;
		temp.vel[Y] = pom*sin(alfa)*sintheta;
		temp.vel[Z] = pom*costheta;
		
		obj->AddStar(gmass,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenGaussStarBall(cObjects *obj, int N, float gmass, float R, float max_vel) {  // koule s normálním rozdělením
	int i;
	float pom,alfa,beta;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		//do {
			temp.pos[X] = R*GenGaussian();
			temp.pos[Y] = R*GenGaussian();
			temp.pos[Z] = R*GenGaussian();
		//} while (sqrt(Rad2(temp.pos))>R);
		
		//do {
			temp.vel[X] = max_vel*GenGaussian();
			temp.vel[Y] = max_vel*GenGaussian();
			temp.vel[Z] = max_vel*GenGaussian();
		//} while (sqrt(Rad2(temp.vel))>max_vel);
		
		obj->AddStar(gmass,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

//--------------------------------------------------------------------

void GenPlummerStar(cObjects *obj, int N, float Plum_gmass, float charR) {  // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
	int i;
	float rad,vel,phi,sintheta,costheta;
	float a,b;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		rad = RandPlummerRad(charR);
		
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.pos[X] = rad*cos(phi)*sintheta;
		temp.pos[Y] = rad*sin(phi)*sintheta;
		temp.pos[Z] = rad*costheta;
		
		vel = RandPlummerVel(Plum_gmass,rad,charR);
				
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.vel[X] = vel*cos(phi)*sintheta;
		temp.vel[Y] = vel*sin(phi)*sintheta;
		temp.vel[Z] = vel*costheta;
		
		obj->AddStar(Plum_gmass/N,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenPlummerSticky(cObjects *obj, int N, float Plum_gmass, float charR, float sticky_gmass, float stickyR) {  // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
	int i;
	float rad,vel,phi,sintheta,costheta;
	float a,b;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		rad = RandPlummerRad(charR);
		
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.pos[X] = rad*cos(phi)*sintheta;
		temp.pos[Y] = rad*sin(phi)*sintheta;
		temp.pos[Z] = rad*costheta;
		
		vel = RandPlummerVel(Plum_gmass,rad,charR);
				
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.vel[X] = vel*cos(phi)*sintheta;
		temp.vel[Y] = vel*sin(phi)*sintheta;
		temp.vel[Z] = vel*costheta;
		
		obj->AddSticky(sticky_gmass,stickyR,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenPlummerStarTrim(cObjects *obj, int N, float Plum_gmass, float charR, float trim) {  // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
	int i;
	float rad,vel,phi,sintheta,costheta;
	float a,b;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		RandPlummerRadVel(Plum_gmass,charR,trim,rad,vel);
		
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.pos[X] = rad*cos(phi)*sintheta;
		temp.pos[Y] = rad*sin(phi)*sintheta;
		temp.pos[Z] = rad*costheta;
		
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.vel[X] = vel*cos(phi)*sintheta;
		temp.vel[Y] = vel*sin(phi)*sintheta;
		temp.vel[Z] = vel*costheta;
		
		obj->AddStar(Plum_gmass/N,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenPlummerStickyTrim(cObjects *obj, int N, float Plum_gmass, float charR, float sticky_gmass, float stickyR, float trim) {  // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
	int i;
	float rad,vel,phi,sintheta,costheta;
	float a,b;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		RandPlummerRadVel(Plum_gmass,charR,trim,rad,vel);
		
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.pos[X] = rad*cos(phi)*sintheta;
		temp.pos[Y] = rad*sin(phi)*sintheta;
		temp.pos[Z] = rad*costheta;
		
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.vel[X] = vel*cos(phi)*sintheta;
		temp.vel[Y] = vel*sin(phi)*sintheta;
		temp.vel[Z] = vel*costheta;
		
		obj->AddSticky(sticky_gmass,stickyR,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenPlummerStarTrim2(cObjects *obj, int N, float Plum_gmass, float charR, float trim) {  // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
	int i;
	float rad,vel,esc,phi,sintheta,costheta;
	float a,b;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		do {
			do rad = RandPlummerRad(charR); while (rad>trim);
			esc = VescPlum(Plum_gmass, rad, charR)-VescPlum(Plum_gmass, trim, charR);
			vel = RandPlummerVel(Plum_gmass,rad,charR);
		} while (vel>esc);
		
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.pos[X] = rad*cos(phi)*sintheta;
		temp.pos[Y] = rad*sin(phi)*sintheta;
		temp.pos[Z] = rad*costheta;
		
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.vel[X] = vel*cos(phi)*sintheta;
		temp.vel[Y] = vel*sin(phi)*sintheta;
		temp.vel[Z] = vel*costheta;
		
		obj->AddStar(Plum_gmass/N,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenPlummerStickyTrim2(cObjects *obj, int N, float Plum_gmass, float charR, float sticky_gmass, float stickyR, float trim) {  // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
	int i;
	float rad,vel,esc,phi,sintheta,costheta;
	float a,b;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		do {
			do rad = RandPlummerRad(charR); while (rad>trim);
			esc = VescPlum(Plum_gmass, rad, charR)-VescPlum(Plum_gmass, trim, charR);
			vel = RandPlummerVel(Plum_gmass,rad,charR);
		} while (vel>esc);
		
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.pos[X] = rad*cos(phi)*sintheta;
		temp.pos[Y] = rad*sin(phi)*sintheta;
		temp.pos[Z] = rad*costheta;
		
		phi = 2*M_PI*GenUniform();
		costheta = 2*GenUniform()-1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.vel[X] = vel*cos(phi)*sintheta;
		temp.vel[Y] = vel*sin(phi)*sintheta;
		temp.vel[Z] = vel*costheta;
		
		obj->AddSticky(sticky_gmass,stickyR,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

//--------------------------------------------------------------------

void GenPlummHommStarDisc(cObjects *obj, int N, float Plum_gmass, float charR, float r, float h) {
	float x, y, z, rad, alfa, vel;
	S_Point temp;
	
	for (int i=0;i<N;i++) {
		rad = r*sqrt(GenUniform());
		alfa = 2*M_PI*GenUniform();
		z = GenUniform(-h,h);
		
		temp.pos[X] = rad*cos(alfa);
		temp.pos[Y] = rad*sin(alfa);
		temp.pos[Z] = z;
		
		vel = sqrt(-rad*AccPlum(Plum_gmass, rad, rad*rad, charR*charR));
		
		temp.vel[X] = -vel*sin(alfa);
		temp.vel[Y] = vel*cos(alfa);
		temp.vel[Z] = 0;
		
		obj->AddStar(Plum_gmass/N,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenPlummHommStickyDisc(cObjects *obj, int N, float Plum_gmass, float charR, float r, float h, float sticky_gmass, float stickyR) {
	float x, y, z, rad, alfa, vel;
	S_Point temp;
	
	for (int i=0;i<N;i++) {
		rad = r*sqrt(GenUniform());
		alfa = 2*M_PI*GenUniform();
		z = GenUniform(-h,h);
		
		temp.pos[X] = rad*cos(alfa);
		temp.pos[Y] = rad*sin(alfa);
		temp.pos[Z] = z;
		
		vel = sqrt(-rad*AccPlum(Plum_gmass, rad, rad*rad, charR*charR));
		
		temp.vel[X] = -vel*sin(alfa);
		temp.vel[Y] = vel*cos(alfa);
		temp.vel[Z] = 0;
		
		obj->AddSticky(sticky_gmass,stickyR,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

//--------------------------------------------------------------------

void GenPlummHommStarDiscHot(cObjects *obj, int N, float Plum_gmass, float charR, float r, float h) {
	float x, y, z, rad, alfa, vel, gauss;
	S_Point temp;
	
	for (int i=0;i<N;i++) {
		rad = r*sqrt(GenUniform());
		alfa = 2*M_PI*GenUniform();
		z = GenUniform(-h,h);
		
		temp.pos[X] = rad*cos(alfa);
		temp.pos[Y] = rad*sin(alfa);
		temp.pos[Z] = z;
		
		gauss = sqrt(-PotPlum(Plum_gmass, rad, charR)*(1-pow(PotPlum(Plum_gmass, r, charR)/PotPlum(Plum_gmass, rad, charR),6))/6);
		
		vel = sqrt(-rad*AccPlum(Plum_gmass, rad, rad*rad, charR*charR)) + gauss*GenGaussian();
		gauss *= GenGaussian();
		
		temp.vel[X] = -vel*sin(alfa)+gauss*cos(alfa);
		temp.vel[Y] = vel*cos(alfa)+gauss*sin(alfa);
		temp.vel[Z] = 0;		
		
		obj->AddStar(Plum_gmass/N,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenPlummHommStickyDiscHot(cObjects *obj, int N, float Plum_gmass, float charR, float r, float h, float sticky_gmass, float stickyR) {
	float x, y, z, rad, alfa, vel, gauss;
	S_Point temp;
	
	for (int i=0;i<N;i++) {
		rad = r*sqrt(GenUniform());
		alfa = 2*M_PI*GenUniform();
		z = GenUniform(-h,h);
		
		temp.pos[X] = rad*cos(alfa);
		temp.pos[Y] = rad*sin(alfa);
		temp.pos[Z] = z;
		
		gauss = sqrt(-PotPlum(Plum_gmass, rad, charR)*(1-pow(PotPlum(Plum_gmass, r, charR)/PotPlum(Plum_gmass, rad, charR),6))/6);
		
		vel = sqrt(-rad*AccPlum(Plum_gmass, rad, rad*rad, charR*charR)) + gauss*GenGaussian();
		gauss *= GenGaussian();
		
		temp.vel[X] = -vel*sin(alfa)+gauss*cos(alfa);
		temp.vel[Y] = vel*cos(alfa)+gauss*sin(alfa);
		temp.vel[Z] = 0;		
		
		obj->AddSticky(sticky_gmass,stickyR,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

//--------------------------------------------------------------------

void GenPlummerStarDisc(cObjects *obj, int N, float Plum_gmass, float charR) {  // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
	int i;
	float rad,vel,phi,sintheta,costheta;
	float a,b;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		rad = RandPlummerRad(charR);
		
		phi = 2*M_PI*GenUniform();
		costheta = 0.2*GenUniform()-0.1;			//hodnoty od -1 do 1 !!! generuje se
		sintheta = sqrt(1-costheta*costheta);	//hodnoty od 0 do 1 !!! počítá se
		
		temp.pos[X] = rad*cos(phi)*sintheta;
		temp.pos[Y] = rad*sin(phi)*sintheta;
		temp.pos[Z] = rad*costheta;
		
		vel = sqrt(-rad*AccPlum(Plum_gmass, rad, rad*rad, charR*charR));
						
		temp.vel[X] = -vel*sin(phi)*sintheta;
		temp.vel[Y] = vel*cos(phi)*sintheta;
		temp.vel[Z] = vel*costheta;
		
		obj->AddStar(Plum_gmass/N,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

void GenPlummerStickyDisc(cObjects *obj, int N, float Plum_gmass, float charR, float stickyR) {  // vygeneruje Plummerovu sféru hvězd s celkovou hmotností Plum_gmass
	int i;
	float rad,vel,phi,sintheta,costheta;
	float a,b;
	S_Point temp;
	
	for (i=0;i<N;i++) {
		rad = RandPlummerRad(charR);
		
		phi = 2*M_PI*GenUniform();
		costheta = 0.2*GenUniform()-0.1;			//hodnoty od -0.1 do 0.1 !!! generuje se výsek z plummera - disk
		sintheta = sqrt(1-costheta*costheta);			//hodnoty od 0 do 1 !!! počítá se
		
		temp.pos[X] = rad*cos(phi)*sintheta;
		temp.pos[Y] = rad*sin(phi)*sintheta;
		temp.pos[Z] = rad*costheta;
		
		vel = sqrt(-rad*AccPlum(Plum_gmass, rad, rad*rad, charR*charR)); //kruhová rychlost + náhodná odchylka
				
		temp.vel[X] = -vel*sin(phi)*sintheta*(0.6*GenUniform()+0.7);
		temp.vel[Y] = vel*cos(phi)*sintheta*(0.6*GenUniform()+0.7);
		temp.vel[Z] = vel*costheta*(0.6*GenUniform()+0.7);
		
		obj->AddSticky(Plum_gmass/N,stickyR,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}



void GenStarLine(cObjects *obj, float gmass, float seg_x, float seg_y, float seg_z, int stars) {
	int i;
	S_Point temp;
	
	temp.pos[X] = 0;	temp.pos[Y] = 0;	temp.pos[Z] = 0;
	temp.vel[X] = 0;	temp.vel[Y] = 0;	temp.vel[Z] = 0;
	
	for (i=0;i<stars;i++) {
		temp.pos[X] += seg_x;
		temp.pos[Y] += seg_y;
		temp.pos[Z] += seg_z;
		
		obj->AddStar(gmass,temp.pos[X],temp.pos[Y],temp.pos[Z],temp.vel[X],temp.vel[Y],temp.vel[Z]);
	}
}

//---------------------------------------------------------------------------------------------------------------------------

void GenDoublePlummer(cObjects *obj, float gmass1, float gmass2, float a1, float a2, int N1, int N2) {

}

void GenDoublePlummer(cObjects *obj, float gmass1, float gmass2, float a1, float a2, int N1, int N2, float trimm1, float trimm2) {

}