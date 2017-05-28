#ifndef _C_DATA_GRID_DEF
	#define _C_DATA_GRID_DEF

#include <vector>
#include <cmath>
#include <fstream>
#include <string>

#define DIM2	2

#define X	0
#define Y	1

class cDataGrid {
	private:
		
	public:
		std::vector<std::vector<float > > data;
		
		float corner1[DIM2];
		float corner2[DIM2];
		
		float size[DIM2];
		float bin_size[DIM2];
		int bins[DIM2];
		
		void SetGrid(float corner1_x, float corner1_y, float corner2_x, float corner2_y, int nbins_x, int nbins_y);
		void AddData(float x, float y, float weight);
		void ClearData();
		
		void SaveGrid(std::string filename);
};

#endif
