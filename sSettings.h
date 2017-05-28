#ifndef _DEF_S_SETTINGS
	#define _DEF_S_SETTINGS

#include <pngwriter.h>

#include <string>

#include "cDataGrid.h"

struct sGraphSetting {
	cDataGrid data_mriz;
	float phys_area[4];
	int bins[2];
	
	pngwriter plot;
	int canv_size[2];
	int edge_size[4];			//rámeček
	std::string png_name;
	
	float tick_dist[2];
	std::string x_axis;
	std::string y_axis;
	
	float ratio[2];		//škálování
	
	float max_value;
	float color[3];
	
	void SetGrid(float min_x, float min_y, float max_x, float max_y, int bins_x, int bins_y);
	void SetPlot(int w, int h, int dx1, int dy1, int dx2, int dy2, std::string name);
	void SetColor(float r, float g, float b, float max);
	void SetAxes(std::string x_name, float x_tick_dist, std::string y_name, float y_tick_dist);
	
	void CalcRatio() {
		ratio[0] = canv_size[0]/(phys_area[2]-phys_area[0]);
		ratio[1] = canv_size[1]/(phys_area[3]-phys_area[1]);
	}
};

struct sGlobalSetting {
	int final_time;
	int draw_every;
	float save;
	int save_every;
	
	int counter;
};

#endif