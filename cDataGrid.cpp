#include "cDataGrid.h"

void cDataGrid::SetGrid(float corner1_x, float corner1_y, float corner2_x, float corner2_y, int nbins_x, int nbins_y) { 
	int i,j;
	
	corner1[X] = corner1_x; corner1[Y] = corner1_y;
	corner2[X] = corner2_x; corner2[Y] = corner2_y;
	
	size[X] = corner2[X] - corner1[X];
	size[Y] = corner2[Y] - corner1[Y];
	bins[X] = nbins_x;
	bins[Y] = nbins_y;
	
	bin_size[X] = size[X]/bins[X];
	bin_size[Y] = size[Y]/bins[Y];
	
	data.resize(bins[X]);
	for (i=0;i<bins[X];i++) data[i].resize(bins[Y]);
	
	ClearData();
}

void cDataGrid::AddData(float x, float y, float weight) {
	int a,b;
	a = trunc((x-corner1[X])/bin_size[X]);
	b = trunc((y-corner1[Y])/bin_size[Y]);
	if (a >= 0 && b >= 0 && a < bins[X] && b < bins[Y]) data[a][b] += weight;
}

void cDataGrid::ClearData() {
	int i,j;
	
	for (i=0;i<bins[X];i++) for (j=0;j<bins[Y];j++) data[i][j] = 0;
}

void cDataGrid::SaveGrid(std::string filename) {
			std::ofstream file;
			file.open(filename.c_str(), std::ios::trunc);
			
			for (int i=0;i<bins[X];i++) for (int j=0;j<bins[Y];j++)
				file << corner1[X]+i*bin_size[X] << "\t" << corner1[Y]+j*bin_size[Y] << "\t" << data[i][j] << std::endl;
		
			file.close();
}