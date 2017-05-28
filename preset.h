#ifndef _C_PRESET
	#define _C_PRESET

#include "params.h"
#include "math_phys.h"
#include "cObjects.h"
#include "cEvolver.h"
#include "sSettings.h"

#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include <pngwriter.h>

#include "GenerateObj.h"
#include "AnalyzeObj.h"
#include "drawing.h"

class cPreset {
	private:
		cObjects *obj;
		cEvolver *evo;
		
		std::ofstream file;
		std::stringstream ss;
		
		std::vector<float > mass_frac;
		std::vector<float > rad;
		std::vector<float > ang_mom;
		std::vector<float > ms_vel;
		
		int pr;
		
		void MergerInit();
		void MergerCtrl();
		void MergerDraw();
		void MergerFinal();
		
		void StickyMergerInit();
		void StickyMergerCtrl();
		void StickyMergerDraw();
		void StickyMergerFinal();
		
		void PlummerInit();
		void PlummerCtrl();
		void PlummerDraw();
		void PlummerFinal();
		
		void StickyPlummerInit();
		void StickyPlummerCtrl();
		void StickyPlummerDraw();
		void StickyPlummerFinal();
		
		void MergerTrimInit();
		void MergerTrimCtrl();
		void MergerTrimDraw();
		void MergerTrimFinal();
		
		void MergerStickyTrimInit();
		void MergerStickyTrimCtrl();
		void MergerStickyTrimDraw();
		void MergerStickyTrimFinal();
		
		void LineInit();
		void LineCtrl();
		void LineDraw();
		void LineFinal();
		
		void HommDiscInit();
		void HommDiscCtrl();
		void HommDiscDraw();
		void HommDiscFinal();
		
		void DiscMergerInit();
		void DiscMergerCtrl();
		void DiscMergerDraw();
		void DiscMergerFinal();
		
		void DynFricTestInit();
		void DynFricTestCtrl();
		void DynFricTestDraw();
		void DynFricTestFinal();
		
		void RozetaInit();
		void RozetaCtrl();
		void RozetaDraw();
		void RozetaFinal();
		
		std::vector<sGraphSetting > gphsett;
		
		void WriteLine(std::string text) {
			file << text << std::flush;
		}
	public:
		sGlobalSetting gsett;
		
		void setting(int pres, cEvolver* e, cObjects* o) {
			obj = o;
			evo = e;
			pr = pres;
		}
		
		void init() {
			switch (pr) {
				case 0:		MergerInit();	break;
				case 1:		StickyMergerInit();	break;
				case 2:		PlummerInit();	break;
				case 3:		StickyPlummerInit();	break;
				case 4:		MergerTrimInit();	break;
				case 5:		MergerStickyTrimInit();	break;
				case 6:		LineInit();	break;
				case 7:		HommDiscInit();	break;
				case 8:		DiscMergerInit();	break;
				case 9:		DynFricTestInit();	break;
				case 10:	RozetaInit();	break;
			}		
		}
		
		void draw() {
			switch (pr) {
				case 0:		MergerDraw();	break;
				case 1:		StickyMergerDraw();	break;
				case 2:		PlummerDraw();	break;
				case 3:		StickyPlummerDraw();	break;
				case 4:		MergerTrimDraw();	break;
				case 5:		MergerStickyTrimDraw();	break;
				case 6:		LineDraw();	break;
				case 7:		HommDiscDraw();	break;
				case 8:		DiscMergerDraw();	break;			
				case 9:		DynFricTestDraw();	break;
				case 10:	RozetaDraw();	break;
			}
		}
		
		void ctrl() {
			switch (pr) {
				case 0:		MergerCtrl();	break;
				case 1:		StickyMergerCtrl();	break;
				case 2:		PlummerCtrl();	break;
				case 3:		StickyPlummerCtrl();	break;
				case 4:		MergerTrimCtrl();	break;
				case 5:		MergerStickyTrimCtrl();	break;
				case 6:		LineCtrl();	break;
				case 7:		HommDiscCtrl();	break;
				case 8:		DiscMergerCtrl();	break;
				case 9:		DynFricTestCtrl();	break;
				case 10:	RozetaCtrl();	break;
			}
		}
		
		void final() {
			switch (pr) {
				case 0:		MergerFinal();	break;
				case 1:		StickyMergerFinal();	break;
				case 2:		PlummerFinal();	break;
				case 3:		StickyPlummerFinal();	break;
				case 4:		MergerTrimFinal();	break;
				case 5:		MergerStickyTrimFinal();	break;
				case 6:		LineFinal();	break;
				case 7:		HommDiscFinal();	break;
				case 8:		DiscMergerFinal();	break;
				case 9:		DynFricTestFinal();	break;
				case 10:	RozetaFinal();	break;
			}
		}
};

#endif