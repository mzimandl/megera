#ifndef _DEF_DRAWFUNC
	#define _DEF_DRAWFUNC

#include <iostream>
#include <sstream>

#include "cObjects.h"
#include "cEvolver.h"
#include "sSettings.h"

#include <pngwriter.h>

#define TICK	2		//délka ticku
#define LABEL_SIZE	6	//velikost popisků
#define TIME_SIZE	12	//velikost textu času
#define X_LABEL_EVERY	2	//popsat každý n-tý x tick
#define Y_LABEL_EVERY	2	//totéž v y
#define X_LABEL_DISP	20	//posunutí popisků
#define Y_LABEL_DISP	20

static const char* font = "unispace.ttf";

void PaletaLog(sGraphSetting &gs, float value, float &r, float &g, float &b);
void PaletaLin(sGraphSetting &gs, float value, float &r, float &g, float &b);
float iPaletaLog(float max, float value);
float iPaletaLin(float max, float value);
void KresliPaletu(sGraphSetting &gs, int x1, int y1, int x2, int y2);
void KresliGrid(sGraphSetting &gs);
void iKresliGrid(sGraphSetting &gs, float r, float g, float b, float mez);
void KresliOsy(sGraphSetting &gs, cObjects *obj);
void UlozPng(sGraphSetting &gs, int index);
void KresliKolizeXY(sGraphSetting &gs, cEvolver *evo);

void KresliPotencial(sGraphSetting &gs, float pos_x, float pos_y, float radius, float opacity, float r, float g, float b);

#endif