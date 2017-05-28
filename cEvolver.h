#ifndef _C_EVOLVER_DEF
	#define _C_EVOLVER_DEF
#include "math_phys.h"
#include <iostream>

#include <vector>
#include <cstddef>
#include <omp.h>
#include "cObjects.h"
#include "units.h"
#include "OMPsettings.h"

struct S_Acc {
	float acc[DIM];
};

struct S_Pos {
	float pos[DIM];
};

class cEvolver{
	private:
		// Modifikované objekty
		cObjects* objects;
		
		// Pomocné proměnné zrychlení
		std::vector<S_Acc > acc_obj;
		std::vector<S_Acc > acc_pot;
		
		std::vector<S_Acc* > acc_thread;
		std::vector<S_StColl* > coll_thread;
		std::vector<int > num_coll_thread;
		
		// Integrační krok		
		float dt;
		
		void ResetAcc();		// Vynuluje zrychlení
		
		bool test_particles, sticky_particles;
		float coll_para, coll_perp, coll_min_vel2;	// Koeficienty kolize sticky - míra snížení vzájemné rychlosti

	public:
		cEvolver();
//		~cEvolver();
		
		int total_num_coll;
		
		void Reset();			// Vynuluje nastavení třídy
		
		void SetTimeStep(float a);	// Nastavení integračního kroku (defaultně 1)
		void SetObjects(cObjects* a);	// Nastavení vyvíjených objektů
		void SetStickyCollParams(float para, float perp, float min_vel);	// Nastavení kolizní koeficienty
		
		void UpdateAcc();		// Výpočet nových gravitačních zrychlení objektů
		void UpdateVel();		// Výpočet nových rychlostí na základě zrychlení
		void UpdatePos();		// Výpočet nových pozic na základě rychlostí
		void HandleStickyCollision();	// Zpracuje kolize sticky particles
		
		//používat pouze společně, od předešlého má jiné indexování pomocných zrychlení
		void UpdateAcc(int G);		// Výpočet nových gravitačních zrychlení objektů dané skupiny
		void UpdateVel(int G);		// Výpočet nových rychlostí na základě zrychlení objektů dané skupiny
		void UpdatePos(int G);		// Výpočet nových pozic na základě rychlostí objektů dané skupiny
		
		void InitLeapFrog();		// Výpočet zrychlení a z něj rychlosti o půl kroku zpět
		void Evolve();			// Update rutiny v pořadí Pos/Acc/Vel + sticky
		void Relax(int G);		// Update rutiny v pořadí Pos/Acc/Vel objetů dané skupiny
		
		void TestParticlesOn() { test_particles = true; };
		void TestParticlesOff() { test_particles = false; };
		
		void StickyParticlesOn() { sticky_particles = true; };
		void StickyParticlesOff() { sticky_particles = false; };
		
		std::vector<S_Pos > coll_pos;
};

#endif
