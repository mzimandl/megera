#include "AnalyzeDataGrid.h"

void FindUpperExtremesRow(cDataGrid *mriz, float x1, float x2, float y, std::vector<float > &output) {
	output.clear();

	x1 = (x1 - mriz->corner1[X])/mriz->bin_size[X];
//	if (x1<0) x1 = 0;
//	else if (x1>mriz->bins[X]) x1 = mriz->bins[X]-1;
	
	x2 = (x2 - mriz->corner1[X])/mriz->bin_size[X];
//	if (x2<0) x2 = 0;
//	else if (x1>mriz->bins[X]) x2 = mriz->bins[X]-1;
	
	y = (y - mriz->corner1[Y])/mriz->bin_size[Y];
//	if (y<0) y = 0;
//	else if (y>mriz->bins[Y]) y = mriz->bins[Y]-1;
	
	int poc = (int)(x1);
	int rozs = (int)(x2-poc);
	int row = (int)(y);
	
	float temp, temp1, temp2;
	
	temp1 = mriz->data[poc+1][row] - mriz->data[poc][row];
	
	for (int i=1;i<rozs-1;i++) {
		//std::cout << temp1 << " | ";
		temp2 = mriz->data[poc+i+1][row] - mriz->data[poc+i][row];
		if (temp1 > 0 && temp2 < 0) {
			temp = mriz->corner1[X]+(poc+i+1/2)*mriz->bin_size[X];
			output.push_back(temp);
		}
		temp1 = temp2;
	} std::cout << std::endl;
	std::cout << "zde" << std::endl;
}


void Convolution(cDataGrid *mriz) {
	float conv_matrix[3][3] = {1,2,1,2,4,2,1,2,1};
	float norm = 0;
	
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) norm += conv_matrix[i][j];
	
	float temp_matrix[mriz->bins[X]][mriz->bins[Y]];
	for (int i=1;i<mriz->bins[X]-1;i++) for (int j=1;j<mriz->bins[Y]-1;j++) {
		temp_matrix[i][j] = 0;
		for (int k=0;k<3;k++) for (int l=0;l<3;l++) temp_matrix[i][j] += mriz->data[i-1+k][j-1+l]*conv_matrix[k][l];
		temp_matrix[i][j] /= norm;
	}
	
	for (int i=1;i<mriz->bins[X]-1;i++) for (int j=1;j<mriz->bins[Y]-1;j++) mriz->data[i][j] = temp_matrix[i][j];
}

void Sobel(cDataGrid *mriz, float prah) {
	
	float conv_matrix[4][3][3] = {-1,-2,-1,0,0,0,1,2,1,-1,0,1,-2,0,2,-1,0,1,-2,-1,0,-1,0,1,0,1,2,0,1,2,-1,0,1,-2,-1,0};
	float temp[4] = {0,0,0,0};
	float max;
	
	float temp_matrix[mriz->bins[X]][mriz->bins[Y]];
	for (int i=1;i<mriz->bins[X]-1;i++) for (int j=1;j<mriz->bins[Y]-1;j++) {
		
		for (int k=0;k<4;k++) temp[k] = 0;
		for (int k=0;k<4;k++) for (int l=0;l<3;l++) for (int m=0;m<3;m++)  temp[k] += mriz->data[i-1+l][j-1+m]*conv_matrix[k][l][m];
		
		max = fabs(temp[0]);
		for (int k=1;k<4;k++) if (fabs(temp[k]) > max) max = fabs(temp[k]);
		if (max > prah) temp_matrix[i][j] = max; else temp_matrix[i][j] = 0;
	}
	
	for (int i=0;i<mriz->bins[X];i++) for (int j=0;j<mriz->bins[Y];j++) mriz->data[i][j] = temp_matrix[i][j];
}

void FastWavelet(cDataGrid &mriz) { //POLE MOCNINY 2!!!
	float zdroj[mriz.bins[X]][mriz.bins[Y]];
	float cil[mriz.bins[X]][mriz.bins[Y]];
	
	for (int i=0;i<mriz.bins[X];i++) for (int j=0;j<mriz.bins[Y];j++) zdroj[i][j] = mriz.data[i][j];
	
	int start, konec;
	
 	start = 0;
	konec = mriz.bins[X]/2;
		
 	while (konec) {
 		for (int i=0;i<konec;i++) for (int j=0;j<mriz.bins[Y];j++) {
			cil[start+i][j] = zdroj[start+2*i][j]-zdroj[start+2*i+1][j];
			cil[start+konec+i][j] = zdroj[start+2*i][j]+zdroj[start+2*i+1][j];
		}
		std::memcpy(zdroj,cil,sizeof(float)*mriz.bins[X]*mriz.bins[Y]);
		
		start += konec;
		konec /= 2;
	}
	
 	start = 0;
	konec = mriz.bins[Y]/2;
		
 	while (konec) {	//opakuje dokud není konec = 0
 		for (int j=0;j<konec;j++) for (int i=0;i<mriz.bins[X];i++) {
			cil[i][start+j] = zdroj[i][start+2*j]-zdroj[i][start+2*j+1];
			cil[i][start+konec+j] = zdroj[i][start+2*j]+zdroj[i][start+2*j+1];
		}
		
		std::memcpy(zdroj,cil,sizeof(float)*mriz.bins[X]*mriz.bins[Y]);
		
		start += konec;
		konec /= 2;
	}
	
	for (int i=0;i<mriz.bins[X];i++) for (int j=0;j<mriz.bins[Y];j++) mriz.data[i][j] = cil[i][j];
}

void ReverseFastWavelet(cDataGrid &mriz) {
	float zdroj[mriz.bins[X]][mriz.bins[Y]];
	float cil[mriz.bins[X]][mriz.bins[Y]];
	
	for (int i=0;i<mriz.bins[X];i++) for (int j=0;j<mriz.bins[Y];j++) zdroj[i][j] = mriz.data[i][j];
	
	int start, konec;
	
 	start = 0;
	konec = mriz.bins[X]/2;
	
 	while (konec) {
 		for (int i=0;i<konec;i++) for (int j=0;j<mriz.bins[Y];j++) {
			cil[start+i][j] = zdroj[start+2*i][j]-zdroj[start+2*i+1][j];
			cil[start+konec+i][j] = zdroj[start+2*i][j]+zdroj[start+2*i+1][j];
		}
		std::memcpy(zdroj,cil,sizeof(float)*mriz.bins[X]*mriz.bins[Y]);
		
		start += konec;
		konec /= 2;
	}
	
 	start = 0;
	konec = mriz.bins[Y]/2;
		
 	while (konec) {
 		for (int j=0;j<konec;j++) for (int i=0;i<mriz.bins[X];i++) {
			cil[i][start+j] = zdroj[i][start+2*j]-zdroj[i][start+2*j+1];
			cil[i][start+konec+j] = zdroj[i][start+2*j]+zdroj[i][start+2*j+1];
		}
		
		std::memcpy(zdroj,cil,sizeof(float)*mriz.bins[X]*mriz.bins[Y]);
		
		start += konec;
		konec /= 2;
	}
	
	for (int i=0;i<mriz.bins[X];i++) for (int j=0;j<mriz.bins[Y];j++) mriz.data[i][j] = cil[i][j];
}