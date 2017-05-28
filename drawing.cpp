#include "drawing.h"

void PaletaLog(sGraphSetting &gs, float value, float &r, float &g, float &b) {
	float intensity = log10(1+9*5*value/gs.max_value);
	if (intensity > 1) intensity = 1;
	else if (intensity < 0) intensity = 0;
	r = intensity*gs.color[0];
	
	intensity = value/gs.max_value;
	if (intensity > 1) intensity = 1;
	else if (intensity < 0) intensity = 0;
	g = intensity*gs.color[1];
	b = intensity*gs.color[2];
}

void PaletaLin(sGraphSetting &gs, float value, float &r, float &g, float &b) {
	float intensity = value/gs.max_value;
	
	if (intensity > 1) intensity = 1;
	else if (intensity < 0) intensity = 0;
	
	r = intensity*gs.color[0];
	g = intensity*gs.color[1];
	b = intensity*gs.color[2];
}

float iPaletaLog(float max, float value) {
	float intensity = log10(1+9*value/max);
	
	if (intensity > 1) intensity = 1;
	else if (intensity < 0) intensity = 0;
	
	return intensity;
}

float iPaletaLin(float max, float value) {
	float intensity = value/max;
	
	if (intensity > 1) intensity = 1;
	else if (intensity < 0) intensity = 0;
	
	return intensity;
}

void KresliPaletu(sGraphSetting &gs, int x1, int y1, int x2, int y2) {
	float delka = y2-y1;
	float temp = gs.max_value/delka;
	float r,g,b;

	for (int i=0;i<delka;i++) {
		PaletaLog(gs, temp*i, r, g, b);
		gs.plot.line(x1,y1+i,x2,y1+i,r,g,b);
	}
	
	gs.plot.square(x1,y1,x2,y2,1.0,1.0,1.0);
	
	std::stringstream ss;
	std::string str;
	
	ss << 0;
	str = ss.str();
	ss.str(std::string());
	
	gs.plot.plot_text((char*)font,LABEL_SIZE,x1,y1-2*LABEL_SIZE,0,(char*)str.c_str(),1.0,1.0,1.0);
	
	ss << gs.max_value;
	str = ss.str();
	ss.str(std::string());
	
	gs.plot.plot_text((char*)font,LABEL_SIZE,x1,y2+LABEL_SIZE,0,(char*)str.c_str(),1.0,1.0,1.0);
}


void KresliGrid(sGraphSetting &gs) {
	float rectx[2], recty[2];
	float r,g,b;
	float intens;
	
	for (int i=0;i<gs.bins[X];i++) for (int j=0;j<gs.bins[Y];j++) {
		rectx[0]=gs.ratio[X]*(gs.data_mriz.corner1[X]+i*gs.data_mriz.bin_size[X]-gs.phys_area[0]);
		rectx[1]=gs.ratio[X]*(gs.data_mriz.corner1[X]+(i+1)*gs.data_mriz.bin_size[X]-gs.phys_area[0]);
		recty[0]=gs.ratio[Y]*(gs.data_mriz.corner1[Y]+j*gs.data_mriz.bin_size[Y]-gs.phys_area[1]);
		recty[1]=gs.ratio[Y]*(gs.data_mriz.corner1[Y]+(j+1)*gs.data_mriz.bin_size[Y]-gs.phys_area[1]);
		
		PaletaLog(gs, gs.data_mriz.data[i][j], r, g, b);
		intens = iPaletaLog(gs.max_value, gs.data_mriz.data[i][j]);
		
		gs.plot.filledsquare_blend(rectx[0]+gs.edge_size[0], recty[0]+gs.edge_size[1], rectx[1]-1+gs.edge_size[0], recty[1]-1+gs.edge_size[1],intens , r, g, b);
	}
}

void iKresliGrid(sGraphSetting &gs, float r, float g, float b, float mez) {
	float rectx[2], recty[2];
	float intens;
	
	for (int i=0;i<gs.bins[X];i++) for (int j=0;j<gs.bins[Y];j++) {
		rectx[0]=gs.ratio[X]*(gs.data_mriz.corner1[X]+i*gs.data_mriz.bin_size[X]-gs.phys_area[0]);
		rectx[1]=gs.ratio[X]*(gs.data_mriz.corner1[X]+(i+1)*gs.data_mriz.bin_size[X]-gs.phys_area[0]);
		recty[0]=gs.ratio[Y]*(gs.data_mriz.corner1[Y]+j*gs.data_mriz.bin_size[Y]-gs.phys_area[1]);
		recty[1]=gs.ratio[Y]*(gs.data_mriz.corner1[Y]+(j+1)*gs.data_mriz.bin_size[Y]-gs.phys_area[1]);
		
		intens = iPaletaLog(mez, gs.data_mriz.data[i][j]);
		
		gs.plot.filledsquare_blend(rectx[0]+gs.edge_size[0], recty[0]+gs.edge_size[1], rectx[1]-1+gs.edge_size[0], recty[1]-1+gs.edge_size[1], intens, r, g, b);
	}
}

void KresliOsy(sGraphSetting &gs, cObjects *obj) {
	int text_w;
	std::stringstream ss;
	std::string str;
	
	gs.plot.square(gs.edge_size[0],gs.edge_size[1],gs.canv_size[0]+gs.edge_size[0],gs.edge_size[1]+gs.canv_size[1],1.0,1.0,1.0);

	gs.plot.line(-gs.ratio[X]*gs.phys_area[0]+gs.edge_size[0],gs.edge_size[1],-gs.ratio[X]*gs.phys_area[0]+gs.edge_size[0],gs.edge_size[1]+4*TICK,1.0,1.0,1.0);
	gs.plot.line(-gs.ratio[X]*gs.phys_area[0]+gs.edge_size[0],gs.edge_size[1]+gs.canv_size[1],-gs.ratio[X]*gs.phys_area[0]+gs.edge_size[0],gs.edge_size[1]-4*TICK+gs.canv_size[1],1.0,1.0,1.0);

	ss << 0;
	str = ss.str();
	ss.str(std::string());
	text_w = gs.plot.get_text_width((char*)font, LABEL_SIZE, (char*)str.c_str());
	gs.plot.plot_text((char*)font,LABEL_SIZE,-gs.ratio[0]*gs.phys_area[0]+gs.edge_size[0]-text_w/2,gs.edge_size[1]-2*LABEL_SIZE,0,(char*)str.c_str(),1.0,1.0,1.0);
	
	for (int i=1;i<=fabs(gs.phys_area[2])/gs.tick_dist[0];i++) {
		gs.plot.line(gs.ratio[X]*(i*gs.tick_dist[0]-gs.phys_area[0])+gs.edge_size[0],gs.edge_size[1],gs.ratio[X]*(i*gs.tick_dist[0]-gs.phys_area[0])+gs.edge_size[0],gs.edge_size[1]+TICK,1.0,1.0,1.0);
		gs.plot.line(gs.ratio[X]*(i*gs.tick_dist[0]-gs.phys_area[0])+gs.edge_size[0],gs.edge_size[1]+gs.canv_size[1],gs.ratio[X]*(i*gs.tick_dist[0]-gs.phys_area[0])+gs.edge_size[0],gs.edge_size[1]-TICK+gs.canv_size[1],1.0,1.0,1.0);
		
		if (i % X_LABEL_EVERY == 0) {
			ss << i*gs.tick_dist[0];
			str = ss.str();
			ss.str(std::string());
//			text_w = gs.plot.get_text_width((char*)font, LABEL_SIZE, (char*)str.c_str());
			gs.plot.plot_text((char*)font,LABEL_SIZE,gs.ratio[X]*(i*gs.tick_dist[0]-gs.phys_area[0])+gs.edge_size[0]-text_w/2,gs.edge_size[1]-2*LABEL_SIZE,0,(char*)str.c_str(),1.0,1.0,1.0);
		}
	}
	
	for (int i=1;i<=fabs(gs.phys_area[0])/gs.tick_dist[0];i++) {
		gs.plot.line(gs.ratio[X]*(-i*gs.tick_dist[0]-gs.phys_area[0])+gs.edge_size[0],gs.edge_size[1],gs.ratio[X]*(-i*gs.tick_dist[0]-gs.phys_area[0])+gs.edge_size[0],gs.edge_size[1]+TICK,1.0,1.0,1.0);
		gs.plot.line(gs.ratio[X]*(-i*gs.tick_dist[0]-gs.phys_area[0])+gs.edge_size[0],gs.edge_size[1]+gs.canv_size[1],gs.ratio[X]*(-i*gs.tick_dist[0]-gs.phys_area[0])+gs.edge_size[0],gs.edge_size[1]-TICK+gs.canv_size[1],1.0,1.0,1.0);
		if (i % X_LABEL_EVERY == 0) {
			ss << -i*gs.tick_dist[0];
			str = ss.str();
			ss.str(std::string());
//			text_w = gs.plot.get_text_width((char*)font, LABEL_SIZE, (char*)str.c_str());
			gs.plot.plot_text((char*)font,LABEL_SIZE,gs.ratio[X]*(-i*gs.tick_dist[0]-gs.phys_area[0])+gs.edge_size[0]-text_w/2,gs.edge_size[1]-2*LABEL_SIZE,0,(char*)str.c_str(),1.0,1.0,1.0);
		}
	}
	
	gs.plot.line(gs.edge_size[0],-gs.ratio[Y]*gs.phys_area[1]+gs.edge_size[1],gs.edge_size[0]+4*TICK,-gs.ratio[Y]*gs.phys_area[1]+gs.edge_size[1],1.0,1.0,1.0);
	gs.plot.line(gs.edge_size[0]+gs.canv_size[0],-gs.ratio[Y]*gs.phys_area[1]+gs.edge_size[1],gs.edge_size[0]+gs.canv_size[0]-4*TICK,-gs.ratio[Y]*gs.phys_area[1]+gs.edge_size[1],1.0,1.0,1.0);
	
	ss << 0;
	str = ss.str();
	ss.str(std::string());
//	text_w = gs.plot.get_text_width((char*)font, LABEL_SIZE, (char*)str.c_str());
	gs.plot.plot_text((char*)font,LABEL_SIZE,gs.edge_size[0]-LABEL_SIZE,-gs.ratio[Y]*gs.phys_area[1]+gs.edge_size[1]-text_w/2,3.14/2,(char*)str.c_str(),1.0,1.0,1.0);
		
	for (int i=1;i<=fabs(gs.phys_area[3])/gs.tick_dist[1];i++) {
		gs.plot.line(gs.edge_size[0],gs.ratio[Y]*(i*gs.tick_dist[1]-gs.phys_area[1])+gs.edge_size[1],gs.edge_size[0]+TICK,gs.ratio[Y]*(i*gs.tick_dist[1]-gs.phys_area[1])+gs.edge_size[1],1.0,1.0,1.0);
		gs.plot.line(gs.edge_size[0]+gs.canv_size[0],gs.ratio[Y]*(i*gs.tick_dist[1]-gs.phys_area[1])+gs.edge_size[1],gs.edge_size[0]+gs.canv_size[0]-TICK,gs.ratio[Y]*(i*gs.tick_dist[1]-gs.phys_area[1])+gs.edge_size[1],1.0,1.0,1.0);
		
		if (i % Y_LABEL_EVERY == 0) {
			ss << i*gs.tick_dist[1];
			str = ss.str();
			ss.str(std::string());
//			text_w = gs.plot.get_text_width((char*)font, LABEL_SIZE, (char*)str.c_str());
			gs.plot.plot_text((char*)font,LABEL_SIZE,gs.edge_size[0]-LABEL_SIZE,gs.ratio[Y]*(i*gs.tick_dist[1]-gs.phys_area[1])+gs.edge_size[1]-text_w/2,3.14/2,(char*)str.c_str(),1.0,1.0,1.0);
		}
	}
	
	for (int i=1;i<=fabs(gs.phys_area[1])/gs.tick_dist[1];i++) {
		gs.plot.line(gs.edge_size[0],gs.ratio[Y]*(-i*gs.tick_dist[1]-gs.phys_area[1])+gs.edge_size[1],gs.edge_size[0]+TICK,gs.ratio[Y]*(-i*gs.tick_dist[1]-gs.phys_area[1])+gs.edge_size[1],1.0,1.0,1.0);
		gs.plot.line(gs.edge_size[0]+gs.canv_size[0],gs.ratio[Y]*(-i*gs.tick_dist[1]-gs.phys_area[1])+gs.edge_size[1],gs.edge_size[0]+gs.canv_size[0]-TICK,gs.ratio[Y]*(-i*gs.tick_dist[1]-gs.phys_area[1])+gs.edge_size[1],1.0,1.0,1.0);
		
		if (i % Y_LABEL_EVERY == 0) {
			ss << -i*gs.tick_dist[1];
			str = ss.str();
			ss.str(std::string());
//			text_w = gs.plot.get_text_width((char*)font, LABEL_SIZE, (char*)str.c_str());
			gs.plot.plot_text((char*)font,LABEL_SIZE,gs.edge_size[0]-LABEL_SIZE,gs.ratio[Y]*(-i*gs.tick_dist[1]-gs.phys_area[1])+gs.edge_size[1]-text_w/2,M_PI/2,(char*)str.c_str(),1.0,1.0,1.0);
		}
	}
			
//	text_w = gs.plot.get_text_width((char*)font, LABEL_SIZE, (char*)gs.x_axis);
	gs.plot.plot_text((char*)font,LABEL_SIZE,gs.canv_size[0]/2,LABEL_SIZE,0,(char*)gs.x_axis.c_str(),1.0,1.0,1.0);
	
//	text_w = gs.plot.get_text_width((char*)font, LABEL_SIZE, (char*)gs.y_axis.c_str());
	gs.plot.plot_text((char*)font,LABEL_SIZE,2*LABEL_SIZE,gs.canv_size[1]/2,M_PI/2,(char*)gs.y_axis.c_str(),1.0,1.0,1.0);
		
	ss << obj->time << " Myr";
	str = ss.str();
	ss.str(std::string());
						
	gs.plot.plot_text((char*)font,TIME_SIZE,10+gs.edge_size[0],10+gs.edge_size[1],0,(char*)str.c_str(),1.0,1.0,1.0);
}

void UlozPng(sGraphSetting &gs, int index) {
	std::stringstream ss;
	std::string str;
	
	ss << index << ".png"; 
	str = gs.png_name + ss.str();
	ss.str(std::string());
	
	gs.plot.pngwriter_rename((char*)str.c_str());
	gs.plot.write_png();
	gs.plot.clear();
}

void KresliKolizeXY(sGraphSetting &gs, cEvolver *evo) {
	for (int i=0; i<evo->coll_pos.size(); i++) if (evo->coll_pos[i].pos[X]>gs.phys_area[0] && evo->coll_pos[i].pos[X]<gs.phys_area[2] && evo->coll_pos[i].pos[Y]>gs.phys_area[1] && evo->coll_pos[i].pos[Y]<gs.phys_area[3])
		gs.plot.plot(gs.edge_size[0]+gs.ratio[X]*(evo->coll_pos[i].pos[X]-gs.phys_area[0]),gs.edge_size[1]+gs.ratio[Y]*(evo->coll_pos[i].pos[Y]-gs.phys_area[1]),0.0,0.0,1.0);
}

void KresliPotencial(sGraphSetting &gs, float pos_x, float pos_y, float radius, float opacity, float r, float g, float b) {
	gs.plot.filledcircle_blend(gs.edge_size[0]+gs.ratio[X]*(pos_x-gs.phys_area[0]), gs.edge_size[1]+gs.ratio[Y]*(pos_y-gs.phys_area[1]), gs.ratio[X]*radius, opacity, r, g, b);
}
