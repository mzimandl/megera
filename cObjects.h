#ifndef _C_OBJECTS_DEF
	#define _C_OBJECTS_DEF

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "math_phys.h"
#include "units.h"

#define O_STAR 0
#define O_STICKY 1

struct S_Group {
   int LStars;
   int LSticky;
   int NStars;
   int NSticky;

   int NStot;
   
   int LPots;
   int NPots;
   
   float relax_time;
};

class cObjects{
	private:
			
	public:
		std::vector<S_Group > group;
		
		//gravitating objects
		std::vector<S_Point > star;
		std::vector<S_Ball > sticky;
		
		int sg_size;				//number of objects
		S_Point* sg_object(int N);		//returns Nth star/sticky object pointer - star structure
		S_Point* sg_object(int G, int N);	//returns Nth star/sticky object pointer from group G - star structure
		
		//solid potentials
		std::vector<S_Potential > potential;
		
		float time;
		
		cObjects();
//		~cObjects();
		
		//resetuje nastavení  
		void Reset();
		
		//přidá skupinu objektů - počet a pozice prvního v poli, první vytvořen automaticky
		void AddGroup();
		
		//přidávají do poslední skupiny příslušné objekty
		void AddStar(float gmass, float x, float y, float z, float vx, float vy, float vz);
		void AddSticky(float gmass, float R , float x, float y, float z, float vx, float vy, float vz);
		void AddPotential(int type, float gmass, float chr, float x, float y, float z, float vx, float vy, float vz);
		
		//úpně vypne vybraný/poslední potenciál
		void TurnOffPotential(int N);
		void TurnOnPotential(int N);
		void TurnOffPotential();
		void TurnOnPotential();
		
		//nastavení počáteku nového souřadného systému
		void SetZeroPoint(float x, float y, float z);
		
		//přičtení souřadnic a rychlostí vybrané/poslední skupině objektů
		void SetGroupOffset(int N, float x, float y, float z, float vx, float vy, float vz);
		void SetGroupOffset(float x, float y, float z, float vx, float vy, float vz);
		
		//Rotace skupiny kolem 0
		void SetGroupRot(int N, float angle, int axis);
		void SetGroupRot(float angle, int axis);
		
		void SetPotentialVel(int N, float vx, float vy, float vz);
		
		void SaveData(std::string filename) {
			std::ofstream file;
			
			std::stringstream ss;
			ss << time;
			filename += ss.str();
			
			file.open(filename.c_str(), std::ios::trunc);
			
			file << "# time " << time << " Myr" << std::endl;
			file << "# groups " << group.size() << " pot/star/sticky" << std::endl;
			for (int i=0;i<group.size();i++) 
				file << "# " << i << "\t" << "relaxed for " << group[i].relax_time << " Myr" << "\t" << group[i].NPots << "\t" << group[i].NStars << "\t" << group[i].NSticky << std::endl;
			
			file << "#" << std::endl;
			file << "# N\tgroup\ttype\tmass[Msun]\tx[kpc]\t\ty[kpc]\t\tz[kpc]\t\tvx[km/s]\tvy[km/s]\tvz[km/s]" << std::endl;
			
			file << std::scientific;
			
			file << "# type 0 - potentials" << std::endl;
			for (int i=0;i<group.size();i++) for (int j=group[i].LPots;j<group[i].LPots+group[i].NPots;j++) {
				file << j << "\t" << i << "\t0\t" << potential[j].gmass*GMUNIT << "\t" << potential[j].pos[X] << "\t" << potential[j].pos[Y] << "\t" << potential[j].pos[Z] << "\t" << potential[j].vel[X]*VUNIT << "\t" << potential[j].vel[Y]*VUNIT << "\t" << potential[j].vel[Z]*VUNIT << std::endl;
			}
			
			file << "# type 1 - stars" << std::endl;
			for (int i=0;i<group.size();i++) for (int j=group[i].LStars;j<group[i].LStars+group[i].NStars;j++) {
				file << j << "\t" << i << "\t1\t" << star[j].gmass*GMUNIT << "\t" << star[j].pos[X] << "\t" << star[j].pos[Y] << "\t" << star[j].pos[Z] << "\t" << star[j].vel[X]*VUNIT << "\t" << star[j].vel[Y]*VUNIT << "\t" << star[j].vel[Z]*VUNIT << std::endl;
			}
			
			file << "# type 2 - sticky" << std::endl;
			for (int i=0;i<group.size();i++) for (int j=group[i].LSticky;j<group[i].LSticky+group[i].NSticky;j++) {
				file << j << "\t" << i << "\t2\t" << sticky[j].gmass*GMUNIT << "\t" << sticky[j].pos[X] << "\t" << sticky[j].pos[Y] << "\t" << sticky[j].pos[Z] << "\t" << sticky[j].vel[X]*VUNIT << "\t" << sticky[j].vel[Y]*VUNIT << "\t" << sticky[j].vel[Z]*VUNIT << std::endl;
			}
			
			file.close();
		}
};

#endif
