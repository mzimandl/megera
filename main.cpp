#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>

#include "cObjects.h"
#include "cEvolver.h"
#include "preset.h"

#define PRESET	9

int main( int argc, const char ** argv ) {
	
	time_t hodiny;		//sytémové hodiny
	hodiny = time(NULL);	//okamžitý čas
	
	srand(hodiny);		//semínko generátoru náhodných čísel

//-------------- deklarace tříd / nastavení výstupu, situace, vývoje -------------------------------  
	
	cPreset preset;		//třída generující specifická nastavení
	
	cObjects objekty;	//informace o částicích
	cEvolver vyvojar;	//integrátor pro třídu cObjects
	
	//seznam nastavení
	preset.setting(PRESET, &vyvojar, &objekty);
	preset.init();
	
	std::cout << 0 << "% " << std::flush;
	
//---------- HLAVNÍ CYKLUS ---------------
	
	int pocet = 0, progress = 1, temp;
	float x, y;
	
	while (objekty.time < preset.gsett.final_time) {
		
		temp = (int)(100*objekty.time/preset.gsett.final_time);
		if (progress != temp) for (progress;progress<=temp;progress++) {
			if (progress%25==0) std::cout << " " << progress << "% " << std::flush;
			else std::cout << "|" << std::flush;
		}
		
		if (pocet % preset.gsett.draw_every == 0) {
		
			preset.draw();
			preset.gsett.counter++;
		}
		
		if (preset.gsett.save && pocet % preset.gsett.save_every == 0) {
			objekty.SaveData("dat");
		}
		
		preset.ctrl();		
		vyvojar.Evolve();	// ZDE SE DĚJE VÝVOJ
	
		pocet++;
	}
	
	std::cout << " 100%\n" << std::flush;

//---------- KONEC HLAVNÍHO CYKLU ---------------
	
	std::stringstream ss;
	std::ofstream myfile;
	myfile.open("cas.txt", std::ios::trunc);
	
	double t = difftime(time(NULL),hodiny);
	ss << "Doba výpočtu: " << (int)(t/3600) << "h ";
	t = t - ((int)(t/3600))*3600;
	ss << (int)(t/60) << "m ";
	t = t - ((int)(t/60))*60;
	ss << t << "s\n";
	
	std::cout << ss.str();
	myfile << ss.str();

	myfile.close();
	
	vyvojar.Reset();
	preset.final();
	
	return 0;
}
