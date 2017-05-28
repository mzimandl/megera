#ifndef _DEF_ANALYZEDATAGRID
	#define _DEF_ANALYZEDATAGRID

#include "cDataGrid.h"
#include <vector>
#include <iostream>
#include <cstring>
#include <cmath>

void FindUpperExtremesRow(cDataGrid *mriz, float x1, float x2, float y, std::vector<float > &output);
void Convolution(cDataGrid *mriz);
void Sobel(cDataGrid *mriz, float prah);
void FastWavelet(cDataGrid &mriz); //POLE MOCNINY 2!!!
void ReverseFastWavelet(cDataGrid &mriz);

#endif