#include "sSettings.h"

void sGraphSetting::SetGrid(float min_x, float min_y, float max_x,  float max_y, int bins_x, int bins_y) {
	phys_area[0] = min_x;
	phys_area[1] = min_y;
	phys_area[2] = max_x;
	phys_area[3] = max_y;
	
	bins[X] = bins_x;
	bins[Y] = bins_y;
	
	data_mriz.SetGrid(min_x, min_y, max_x, max_y, bins_x, bins_y);
};

void sGraphSetting::SetPlot(int w, int h, int dx1, int dy1, int dx2, int dy2, std::string name) {
	canv_size[X] = w;
	canv_size[Y] = h;
	
	edge_size[0] = dx1;
	edge_size[1] = dy1;
	edge_size[2] = dx2;
	edge_size[3] = dy2;
	
	png_name = name;
	plot = pngwriter(w+dx1+dx2,h+dy1+dy2,0.0,"");
};

void sGraphSetting::SetColor(float r, float g, float b, float max) {
	color[0] = r;
	color[1] = g;
	color[2] = b;
	
	max_value = max;
};

void sGraphSetting::SetAxes(std::string x_name, float x_tick_dist, std::string y_name, float y_tick_dist) {
	tick_dist[X] = x_tick_dist;
	x_axis = x_name;
	
	tick_dist[Y] = y_tick_dist;
	y_axis = y_name;
};
