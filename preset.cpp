#include "preset.h"

void cPreset::MergerInit() {
	file.open("analyza.txt", std::ios::trunc);
	
	gsett.final_time = 2001;
	gsett.draw_every = 2000;
	gsett.save_every = 2000;
	gsett.save = true;
	
	gsett.counter = 0;
	
	gphsett.resize(3);
	
	gphsett[0].SetGrid(-50, -50, 50, 50, 400, 400);
	gphsett[0].SetPlot(400, 400, 40, 40, 20, 20, "dens");
	gphsett[0].SetAxes("x [kpc]", 10, "y [kpc]", 10);
	gphsett[0].SetColor(1.0,1.0,0.0,150);
	gphsett[0].CalcRatio();
		
	gphsett[1].SetGrid(-50, -800, 50, 800, 400, 400);
	gphsett[1].SetPlot(400, 400, 40, 40, 20, 20, "x_vx");
	gphsett[1].SetAxes("x [kpc]", 10, "vx [km/s]", 100);
	gphsett[1].SetColor(1.0,1.0,0.0,150);
	gphsett[1].CalcRatio();
	
	gphsett[2].SetGrid(-50, -400, 50, 400, 400, 400);
	gphsett[2].SetPlot(400, 400, 40, 40, 20, 20, "x_vz");
	gphsett[2].SetAxes("x [kpc]", 10, "vz [km/s]", 100);
	gphsett[2].SetColor(1.0,1.0,0.0,150);
	gphsett[2].CalcRatio();

	ss.str(std::string());
	float escape_vel = sqrt(-2*PotPlum(GMP, Dini, PPlumSc)-2*PotPlum(GMS, Dini, SPlumSc));
	std::cout << escape_vel*VUNIT << "\n" << GMStoP*escape_vel*VUNIT/(1+GMStoP) << "\n" << escape_vel*VUNIT/(1+GMStoP) << "\n";
	
	ss.str(std::string());
	ss << "#plumer1 " << GMP*GMUNIT << " Msun, scale " << PPlumSc << " kpc" << std::endl;
	ss << "#plumer2 " << GMS*GMUNIT << " Msun, scale " << SPlumSc << " kpc" << std::endl;
	ss << "#3000000 stars" << std::endl;
	ss << "#esc vel: " << escape_vel*VUNIT << ", sec " << GMStoP*escape_vel*VUNIT/(1+GMStoP) << ", prim " << escape_vel*VUNIT/(1+GMStoP) << " (km/s)" << std::endl;
	
	const float rel_vel = 0.75;			// > 1 - unbound, = 1 - parabolic, < 1 - bound
	const float rel_x_vel = 1.0;
				
	escape_vel *= rel_vel;
	obj->AddPotential(T_PLUM,GMP,PPlumSc,0,0,0,0,0,0);
	GenPlummerStar(obj,0,GMP,PPlumSc);
	obj->SetGroupOffset(0,0,0,-GMStoP*rel_x_vel*escape_vel/(1+GMStoP),-GMStoP*sqrt(1-rel_x_vel*rel_x_vel)*escape_vel/(1+GMStoP),0);
	
	obj->AddGroup();
	obj->AddPotential(T_PLUM,GMS,SPlumSc,0,0,0,0,0,0);
	GenPlummerStar(obj,1000000,GMS,SPlumSc);
	obj->SetGroupOffset(-DTiniS-DTiniP,0,0,rel_x_vel*escape_vel/(1+GMStoP),sqrt(1-rel_x_vel*rel_x_vel)*escape_vel/(1+GMStoP),0);
	
	ss << "#rel vel: " << escape_vel*VUNIT << ", sec " << GMStoP*escape_vel*VUNIT/(1+GMStoP) << ", prim " << escape_vel*VUNIT/(1+GMStoP) << " (km/s)" << std::endl;
	ss << "#separation: " << Dini << " kpc" << std::endl;
//	ss << "#T [Myr]\tN_coll" << std::endl;
	WriteLine(ss.str());
	
	evo->SetObjects(obj);
	evo->InitLeapFrog();
}

void cPreset::MergerDraw() {
	obj->SetZeroPoint(obj->potential[0].pos[X],obj->potential[0].pos[Y],0);
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], 0, obj->star.size());
	FillPhaseGrid(obj, gphsett[1], 0, obj->star.size(),0,0);
	FillPhaseGrid(obj, gphsett[2], 0, obj->star.size(),0,2);
	
	for (int i=0;i<gphsett.size();i++) {
		KresliGrid(gphsett[i]);
		KresliOsy(gphsett[i], obj);
		KresliPaletu(gphsett[i],gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+2, gphsett[i].edge_size[1], gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+gphsett[i].edge_size[2]-2, gphsett[i].canv_size[Y]+gphsett[i].edge_size[1]);
		UlozPng(gphsett[i], gsett.counter);
	}
	
	if (obj->time == 2000) gphsett[0].data_mriz.SaveGrid("grid2000");

}

void cPreset::MergerCtrl() {
	if (obj->potential[1].pos[X] >= obj->potential[0].pos[X] && obj->potential[1].active) {
		obj->SetPotentialVel(0,0,0,0);
		obj->SetPotentialVel(1,0,0,0);
		
		obj->TurnOffPotential(1);
	}
	
//	ss.str(std::string());
//	ss << obj->time << std::endl;
//	WriteLine(ss.str());
}

void cPreset::MergerFinal() {
	file.close();
	for (int i=0;i<gphsett.size();i++) gphsett[i].plot.close();
}

//--------------------------------------------------------------------------------------------------

void cPreset::StickyMergerInit() {
	const float rel_vel = 0.75;			// > 1 - unbound, = 1 - parabolic, < 1 - bound
	const float rel_x_vel = 1.0;
	const float separace = 100;
	
	file.open("analyza.txt", std::ios::trunc);
	
	gsett.final_time = 1001;
	gsett.draw_every = 2;
	gsett.save = true;
	gsett.save_every = 100;
	
	gsett.counter = 0;
	
	gphsett.resize(4);
	
	gphsett[0].SetGrid(-50, -50, 50, 50, 200, 200);
	gphsett[0].SetPlot(400, 400, 40, 40, 20, 20, "dens");
	gphsett[0].SetAxes("x [kpc]", 10, "y [kpc]", 10);
	gphsett[0].SetColor(1.0,1.0,0.0,30);
	gphsett[0].CalcRatio();
		
	gphsett[1].SetGrid(-50, -800, 50, 800, 200, 200);
	gphsett[1].SetPlot(400, 400, 40, 40, 20, 20, "x_vx");
	gphsett[1].SetAxes("x [kpc]", 10, "vx [km/s]", 100);
	gphsett[1].SetColor(1.0,1.0,0.0,30);
	gphsett[1].CalcRatio();
	
	gphsett[2].SetGrid(-50, -400, 50, 400, 200, 200);
	gphsett[2].SetPlot(400, 400, 40, 40, 20, 20, "x_vz");
	gphsett[2].SetAxes("x [kpc]", 10, "vz [km/s]", 100);
	gphsett[2].SetColor(1.0,1.0,0.0,30);
	gphsett[2].CalcRatio();
	
	gphsett[3].SetGrid(-100, -50, 0, 50, 200, 200);
	gphsett[3].SetPlot(400, 400, 40, 40, 20, 20, "small");
	gphsett[3].SetAxes("x [kpc]", 10, "y [kpc]", 10);
	gphsett[3].SetColor(1.0,1.0,0.0,30);
	gphsett[3].CalcRatio();
	
	float escape_vel = sqrt(-2*PotPlum(GMP, separace, PPlumSc)-2*PotPlum(GMS, separace, SPlumSc));
	std::cout << escape_vel*VUNIT << "\n" << GMStoP*escape_vel*VUNIT/(1+GMStoP) << "\n" << escape_vel*VUNIT/(1+GMStoP) << "\n";
	
	ss.str(std::string());
	ss << "#plumer1 " << GMP*GMUNIT << " Msun, scale " << PPlumSc << " kpc" << std::endl;
	ss << "#plumer2 " << GMS*GMUNIT << " Msun, scale " << SPlumSc << " kpc" << std::endl;
	ss << "#100000 stars, 10000 sticky - r: 0.05 kpc, Cpara 0.5, Cperp 1.0" << std::endl;
	ss << "#esc vel " << escape_vel*VUNIT << ", sec " << GMStoP*escape_vel*VUNIT/(1+GMStoP) << ", prim " << escape_vel*VUNIT/(1+GMStoP) << " (km/s)" << std::endl;
				
	escape_vel *= rel_vel;
	obj->AddPotential(T_PLUM,GMP,PPlumSc,0,0,0,0,0,0);
	GenPlummerStar(obj,0,GMP,PPlumSc);
	obj->SetGroupOffset(0,0,0,-GMStoP*rel_x_vel*escape_vel/(1+GMStoP),-GMStoP*sqrt(1-rel_x_vel*rel_x_vel)*escape_vel/(1+GMStoP),0);
	
	obj->AddGroup();
	obj->AddPotential(T_PLUM,GMS,SPlumSc,0,0,0,0,0,0);
	GenPlummerStar(obj,100000,GMS,SPlumSc);
	GenPlummerSticky(obj,10000,GMS,SPlumSc,1000/GMUNIT,0.05);
	obj->SetGroupOffset(-separace,0,0,rel_x_vel*escape_vel/(1+GMStoP),sqrt(1-rel_x_vel*rel_x_vel)*escape_vel/(1+GMStoP),0);
	
	evo->SetStickyCollParams(0.5,1.0,1);
	
	ss << "#rel vel " << escape_vel*VUNIT << ", sec " << GMStoP*escape_vel*VUNIT/(1+GMStoP) << ", prim " << escape_vel*VUNIT/(1+GMStoP) << " (km/s)" << std::endl;
	ss << "#separation " << Dini << " kpc" << std::endl;
	ss << "#GMUNIT " << GMUNIT << ", VUNIT" << VUNIT << std::endl;
	ss << "#T [Myr]\tN_coll\tAngM/M 0.1kpc\tAngM/M 0.5kpc\tAngM/M 1kpc\tAngM/M 2kpc\tAngM/M 5kpc\tAngM/M 10kpc\tAngM/M 50kpc\tN 0.1kpc\tN 0.5kpc\tN 1kpc\tN 2kpc\tN 5kpc\tN 10kpc\tN 50kpc" << std::endl;
	WriteLine(ss.str());
	
	evo->SetObjects(obj);
	evo->InitLeapFrog();
	evo->StickyParticlesOn();
}

void cPreset::StickyMergerDraw() {
	obj->SetZeroPoint(obj->potential[0].pos[X],obj->potential[0].pos[Y],0);
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], 0, obj->star.size());
	FillPhaseGrid(obj, gphsett[1], 0, obj->star.size(),0,0);
	FillPhaseGrid(obj, gphsett[2], 0, obj->star.size(),0,2);
	FillDensityGrid(obj, gphsett[3], 0, obj->star.size());
	
	for (int i=0;i<gphsett.size();i++) {
		KresliGrid(gphsett[i]);
	}	
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], obj->star.size(), obj->star.size()+obj->sticky.size());
	FillPhaseGrid(obj, gphsett[1], obj->star.size(), obj->star.size()+obj->sticky.size(),0,0);
	FillPhaseGrid(obj, gphsett[2], obj->star.size(), obj->star.size()+obj->sticky.size(),0,2);	
	FillDensityGrid(obj, gphsett[3], obj->star.size(), obj->star.size()+obj->sticky.size());
		
	for (int i=0;i<gphsett.size();i++) {
		iKresliGrid(gphsett[i],0.0,0.0,1.0,15);
		KresliOsy(gphsett[i], obj);
//		KresliPaletu(gphsett[i],gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+2, gphsett[i].edge_size[1], gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+gphsett[i].edge_size[2]-2, gphsett[i].canv_size[Y]+gphsett[i].edge_size[1]);
		UlozPng(gphsett[i], gsett.counter);
	}

}

void cPreset::StickyMergerCtrl() {
	ss.str(std::string());
	ss << obj->time << "\t" << evo->total_num_coll << "\t";
	ss << AngMomentumMag(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 0.1) << "\t";
	ss << AngMomentumMag(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 0.5) << "\t";
	ss << AngMomentumMag(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 1) << "\t";
	ss << AngMomentumMag(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 2) << "\t";
	ss << AngMomentumMag(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 5) << "\t";
	ss << AngMomentumMag(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 10) << "\t";
	ss << AngMomentumMag(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 50) << "\t";
	
	ss << NParticles(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 0.1) << "\t";
	ss << NParticles(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 0.5) << "\t";
	ss << NParticles(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 1) << "\t";
	ss << NParticles(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 2) << "\t";
	ss << NParticles(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 5) << "\t";
	ss << NParticles(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 10) << "\t";
	ss << NParticles(obj,obj->star.size(), obj->star.size()+obj->sticky.size(), &obj->potential[1], 50) << std::endl;
	
	
	if (obj->potential[1].active) {
			if (obj->potential[1].pos[X] >= obj->potential[0].pos[X] && obj->potential[1].active) {
			obj->SetPotentialVel(0,0,0,0);
			obj->SetPotentialVel(1,0,0,0);
		
			obj->TurnOffPotential(1);
		}
	}
	
	WriteLine(ss.str());
}

void cPreset::StickyMergerFinal() {
	file.close();
	for (int i=0;i<gphsett.size();i++) gphsett[i].plot.close();
}

//--------------------------------------------------------------------------------------------------

void cPreset::PlummerInit() {
	file.open("analyza.txt", std::ios::trunc);
	
	gsett.final_time = 1001;
	gsett.draw_every = 1;
	gsett.save = false;
	gsett.counter = 0;
	
	gphsett.resize(1);
	gphsett[0].SetGrid(-2, -2, 2, 2, 200, 200);
	gphsett[0].SetPlot(201, 201, 40, 40, 20, 20, "dens");
	gphsett[0].SetAxes("x [kpc]", 0.5, "y [kpc]", 0.5);
	gphsett[0].SetColor(1.0,1.0,0.0,250);
	gphsett[0].CalcRatio();
		
	for(int i=1; i<=10; i++) mass_frac.push_back(i*1.0/10);
	rad.resize(mass_frac.size());
	ang_mom.resize(mass_frac.size());
	ms_vel.resize(mass_frac.size());
	
	ss.str(std::string());
	ss << "#Plummer " << GMS*GMUNIT << " M_sun, scale " << SPlumSc << " kpc" << std::endl;
	ss << "#1000000 stars, trimm2 " << 10*SPlumSc  << " plummer, scale " << SPlumSc  << std::endl;
	ss << "#GMUNIT [M_sun] " << GMUNIT << ", VUNIT [km/s]" << VUNIT << std::endl;
	
	ss << "#T [Myr]\tN_coll";
	for (int i=1; i<=mass_frac.size(); i++) ss <<"\t|" << i << "\tLagM\tLagR [kpc]\tang_mom [kpc^2/Myr]\tmean_square_vel [kpc^2/Myr^2]";
	ss << std::endl;
	WriteLine(ss.str());
	
	obj->AddPotential(T_PLUM,GMS,SPlumSc,0,0,0,0,0,0);
	GenPlummerStarTrim2(obj,1000000,GMS,SPlumSc,10*SPlumSc);
	//GenPlummerStar(obj,1000000,GMS,SPlumSc);
	
	evo->SetObjects(obj);
	evo->InitLeapFrog();
	evo->StickyParticlesOff();
}

void cPreset::PlummerDraw() {
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], 0, obj->star.size());
		
	for (int i=0;i<gphsett.size();i++) {
		KresliGrid(gphsett[i]);
		KresliOsy(gphsett[i], obj);
//		KresliPaletu(gphsett[i],gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+2, gphsett[i].edge_size[1], gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+gphsett[i].edge_size[2]-2, gphsett[i].canv_size[Y]+gphsett[i].edge_size[1]);
		UlozPng(gphsett[i], gsett.counter);
	}

}

void cPreset::PlummerCtrl() {
	
	StarLagrangeR(obj,&obj->potential[0],mass_frac,rad);
	StarAngMomentumMag(obj,&obj->potential[0],rad,ang_mom);
	StarMeanSquareVel2Mag(obj,&obj->potential[0],rad,ms_vel);
	
	ss.str(std::string());
	ss << obj->time << "\t" << evo->total_num_coll;
	for (int i=0; i<mass_frac.size(); i++) ss <<"\t|\t" << mass_frac[i] << "\t" << rad[i] << "\t" << ang_mom[i] << "\t" << ms_vel[i];
	ss << std::endl;
	
	WriteLine(ss.str());
}

void cPreset::PlummerFinal() {
	file.close();
	for (int i=0;i<gphsett.size();i++) gphsett[i].plot.close();
}

//--------------------------------------------------------------------------------------------------

void cPreset::StickyPlummerInit() {
	file.open("analyza.txt", std::ios::trunc);
	
	gsett.final_time = 2001;
	gsett.draw_every = 2000;
	gsett.save = true;
	gsett.save_every = 500;
	
	gsett.counter = 0;
	
	gphsett.resize(1);
	
	gphsett[0].SetGrid(-2.5, -2.5, 2.5, 2.5, 400, 400);
	gphsett[0].SetPlot(400, 400, 40, 40, 20, 20, "dens");
	gphsett[0].SetAxes("x [kpc]", 1.0, "y [kpc]", 1.0);
	gphsett[0].SetColor(1.0,1.0,0.0,5);
	gphsett[0].CalcRatio();
		
/*	gphsett[1].SetGrid(-100, -100, 100, 100, 200, 200);
	gphsett[1].SetPlot(400, 400, 40, 40, 20, 20, "vx_vy");
	gphsett[1].SetAxes("vx [km/s]", 20, "vy [km/s]", 20);
	gphsett[1].SetColor(1.0,1.0,0.0,15);
	gphsett[1].CalcRatio();
*/	
	for(int i=1; i<=10; i++) mass_frac.push_back(i*1.0/10);
	rad.resize(mass_frac.size());
	ang_mom.resize(mass_frac.size());
	ms_vel.resize(mass_frac.size());
	
	ss.str(std::string());
	ss << "#plumer " << GMS*GMUNIT << " Msun, scale " << SPlumSc << " kpc, trimm" << 10*SPlumSc << " kpc" << std::endl;
	ss << "#10000 sticky - r: 0.020 kpc, Cpara 0.5, Cperp 1.0 " << std::endl;
	
	obj->AddPotential(T_PLUM,GMS,SPlumSc,0,0,0,0,0,0);
//	GenPlummerSticky(obj,20000,GMS,SPlumSc,2000/GMUNIT,0.010);
	GenPlummerStickyTrim2(obj,10000,GMS,SPlumSc,2000/GMUNIT,0.015,10*SPlumSc);
	
	evo->SetStickyCollParams(0.5,1.0,5);
	
	ss << "#GMUNIT " << GMUNIT << ", VUNIT" << VUNIT << std::endl;
	ss << "#T [Myr]\tN_coll";
	for (int i=1; i<=mass_frac.size(); i++) ss <<"\t|" << i << "\tLagM\tLagR [kpc]\tang_mom [kpc^2/Myr]\tmean_square_vel [kpc^2/Myr^2]";
	ss << std::endl;
	WriteLine(ss.str());
	
	evo->SetObjects(obj);
	evo->InitLeapFrog();
	evo->StickyParticlesOff();
	
	std::cout << "relaxuji..." << std::endl;
	for (int i=0; i<1000; i++) evo->Evolve();
	
	obj->time = 0;
	evo->StickyParticlesOn();
}

void cPreset::StickyPlummerDraw() {
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], obj->star.size(), obj->star.size()+obj->sticky.size());
//	FillVelocityGrid(obj, gphsett[1], obj->star.size(), obj->star.size()+obj->sticky.size());
		
	for (int i=0;i<gphsett.size();i++) {
		KresliGrid(gphsett[i]);
		KresliOsy(gphsett[i], obj);
//		KresliPaletu(gphsett[i],gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+2, gphsett[i].edge_size[1], gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+gphsett[i].edge_size[2]-2, gphsett[i].canv_size[Y]+gphsett[i].edge_size[1]);
		UlozPng(gphsett[i], gsett.counter);
	}

}

void cPreset::StickyPlummerCtrl() {
	StickyLagrangeR(obj,&obj->potential[0],mass_frac,rad);
	StickyAngMomentumMag(obj,&obj->potential[0],rad,ang_mom);
	StickyMeanSquareVel2Mag(obj,&obj->potential[0],rad,ms_vel);
	
	ss.str(std::string());
	ss << obj->time << "\t" << evo->total_num_coll;
	for (int i=0; i<mass_frac.size(); i++) ss <<"\t|\t" << mass_frac[i] << "\t" << rad[i] << "\t" << ang_mom[i] << "\t" << ms_vel[i];
	ss << std::endl;
	
	WriteLine(ss.str());
}

void cPreset::StickyPlummerFinal() {
	file.close();
	for (int i=0;i<gphsett.size();i++) gphsett[i].plot.close();
}

//---------------------------------------------------------------------------------------------------------------------------------

void cPreset::MergerTrimInit() {
	file.open("analyza.txt", std::ios::trunc);
	
	gsett.final_time = 2001;
	gsett.draw_every = 2000;
	gsett.save = true;
	gsett.save_every = 2000;
	
	gsett.counter = 0;
	
	gphsett.resize(1);
	
	gphsett[0].SetGrid(-50, -50, 50, 50, 400, 400);
	gphsett[0].SetPlot(400, 400, 40, 40, 20, 20, "dens");
	gphsett[0].SetAxes("x [kpc]", 10, "y [kpc]", 10);
	gphsett[0].SetColor(0,0,1.0,50);
	gphsett[0].CalcRatio();
/*		
	gphsett[1].SetGrid(-50, -800, 50, 800, 400, 400);
	gphsett[1].SetPlot(400, 400, 40, 40, 20, 20, "x_vx");
	gphsett[1].SetAxes("x [kpc]", 10, "vx [km/s]", 100);
	gphsett[1].SetColor(0,0,1.0,50);
	gphsett[1].CalcRatio();
	
	gphsett[2].SetGrid(-50, -400, 50, 400, 400, 400);
	gphsett[2].SetPlot(400, 400, 40, 40, 20, 20, "x_vz");
	gphsett[2].SetAxes("x [kpc]", 10, "vz [km/s]", 100);
	gphsett[2].SetColor(0,0,1.0,50);
	gphsett[2].CalcRatio();
*/
	ss.str(std::string());
	float escape_vel = sqrt(-2*PotPlum(GMP, Dini, PPlumSc)-2*PotPlum(GMS, Dini, SPlumSc));
	std::cout << escape_vel*VUNIT << "\n" << GMStoP*escape_vel*VUNIT/(1+GMStoP) << "\n" << escape_vel*VUNIT/(1+GMStoP) << "\n";
	
	ss.str(std::string());
	ss << "#plumer1 " << GMP*GMUNIT << " Msun, scale " << PPlumSc << " kpc" << std::endl;
	ss << "#plumer2 " << GMS*GMUNIT << " Msun, scale " << SPlumSc << " kpc" << std::endl;
	ss << "#1000000 stars - SPlumSc, trim 5*SPlumSc; 1000000 stars - SPlumSc, trim 10*SPlumSc" << std::endl;
	ss << "#esc vel: " << escape_vel*VUNIT << ", sec " << GMStoP*escape_vel*VUNIT/(1+GMStoP) << ", prim " << escape_vel*VUNIT/(1+GMStoP) << " (km/s)" << std::endl;
//	ss << "#s dynamickym trenim - half mass primaru - core radius sekundaru, neradialni vx=0.95v";
	
	const float rel_vel = 0.75;			// > 1 - unbound, = 1 - parabolic, < 1 - bound
	const float rel_x_vel = 1.0;
				
	escape_vel *= rel_vel;
	obj->AddPotential(T_PLUM,GMP,PPlumSc,0,0,0,0,0,0);
	//GenPlummerStar(obj,0,GMP,PPlumSc);
	
	obj->AddGroup();
	obj->AddPotential(T_PLUM,GMS,SPlumSc,0,0,0,0,0,0);
	GenPlummerStarTrim2(obj,1000000,GMS,SPlumSc,5*SPlumSc);
	GenPlummerStarTrim2(obj,1000000,GMS,SPlumSc,10*SPlumSc);
	
	
	ss << "#rel vel: " << escape_vel*VUNIT << ", sec " << GMStoP*escape_vel*VUNIT/(1+GMStoP) << ", prim " << escape_vel*VUNIT/(1+GMStoP) << " (km/s)" << std::endl;
	ss << "#separation: " << Dini << " kpc" << std::endl;
//	ss << "#T [Myr]\tDynFric X\tDynFric Y\tDynFric Z" << std::endl;
	WriteLine(ss.str());
	
	evo->SetObjects(obj);
	evo->InitLeapFrog();

/*	
	sGraphSetting relaxgph;
	relaxgph.SetGrid(-10, -10, 10, 10, 400, 400);
	relaxgph.SetPlot(400, 400, 40, 40, 20, 20, "relax");
	relaxgph.SetAxes("x [kpc]", 2, "y [kpc]", 2);
	relaxgph.SetColor(0,0,1.0,40);
	relaxgph.CalcRatio();
*/
	
	int count=0;
	std::cout << "relaxuji..." << std::endl;
	for (int i=0; i<1000; i++) { //relaxace
/*		if (i%2 == 0) {
			relaxgph.data_mriz.ClearData();
			FillDensityGrid(obj, relaxgph, 1000000, obj->star.size());
			iKresliGrid(relaxgph,1.0,1.0,0,40);
			
			relaxgph.data_mriz.ClearData();
			FillDensityGrid(obj, relaxgph, 0, 1000000);
			KresliGrid(relaxgph);
			
			KresliOsy(relaxgph, obj);
			
			UlozPng(relaxgph, count);
			
			count++;
		}
*/
		evo->Relax(1);
	}
	
	obj->time = 0;
	
	obj->SetGroupOffset(0,0,0,0,-GMStoP*rel_x_vel*escape_vel/(1+GMStoP),-GMStoP*sqrt(1-rel_x_vel*rel_x_vel)*escape_vel/(1+GMStoP),0);
	obj->SetGroupOffset(1,-DTiniS-DTiniP,0,0,rel_x_vel*escape_vel/(1+GMStoP),sqrt(1-rel_x_vel*rel_x_vel)*escape_vel/(1+GMStoP),0);
	
	std::cout << "simuluji..." << std::endl;
}

void cPreset::MergerTrimDraw() {
	obj->SetZeroPoint(obj->potential[0].pos[X],obj->potential[0].pos[Y],0);
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], 1000000, obj->star.size());
//	FillPhaseGrid(obj, gphsett[1], 1000000, obj->star.size(),0,0);
//	FillPhaseGrid(obj, gphsett[2], 1000000, obj->star.size(),0,2);	

	for (int i=0;i<gphsett.size();i++) {
		iKresliGrid(gphsett[i],1.0,1.0,0,50);
	}	
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], 0, 1000000);
//	FillPhaseGrid(obj, gphsett[1], 0, 1000000,0,0);
//	FillPhaseGrid(obj, gphsett[2], 0, 1000000,0,2);
			
	for (int i=0;i<gphsett.size();i++) {
		KresliGrid(gphsett[i]);
		KresliOsy(gphsett[i], obj);
//		KresliPaletu(gphsett[i],gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+2, gphsett[i].edge_size[1], gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+gphsett[i].edge_size[2]-2, gphsett[i].canv_size[Y]+gphsett[i].edge_size[1]);
		UlozPng(gphsett[i], gsett.counter);
	}
}

void cPreset::MergerTrimCtrl() {
	
	float radvel=0;
	
	if (obj->potential[1].active) {
		for (int i=0; i<3; i++ ) radvel += (obj->potential[1].vel[i]-obj->potential[0].vel[i])*(obj->potential[1].pos[i]-obj->potential[0].pos[i]);
		
		if (radvel > 0) {
			obj->SetPotentialVel(0,0,0,0);
			obj->SetPotentialVel(1,0,0,0);
		
			obj->TurnOffPotential(1);
			
			ss << "rozpad: " << obj->time << std::endl;
			WriteLine(ss.str());
		}
	}
/*	
	float acc[3];
	
	if (obj->potential[1].active) {
		HandleDynamicalFrictionPlummer(obj->potential[0], obj->potential[1], acc);
	
		for (int i=0;i<3;i++) {
			if (abs(obj->potential[1].vel[i])<abs(acc[i]) && obj->potential[1].vel[i]*acc[i]<0) obj->potential[1].vel[i]=0;
			else obj->potential[1].vel[i]+=acc[i]*1;
//			std::cout << acc[i] << "	";
		}
//		std::cout << std::endl;
	}
*/
/*	
	ss.str(std::string());
	ss << obj->time << "\t" << acc[0] << "\t" << acc[1] << "\t" << acc[2] << std::endl;
	WriteLine(ss.str());
*/
}

void cPreset::MergerTrimFinal() {
	file.close();
	for (int i=0;i<gphsett.size();i++) gphsett[i].plot.close();
}

//------------------------------------------------------------------------------------------------------

void cPreset::MergerStickyTrimInit() {
	const float rel_vel = 0.75;			// > 1 - unbound, = 1 - parabolic, < 1 - bound
	const float rel_x_vel = 1.0;
	const float separace = 100;
	
	file.open("analyza.txt", std::ios::trunc);
	
	gsett.final_time = 2001;
	gsett.draw_every = 2000;
	gsett.save = true;
	gsett.save_every = 2000;
	
	gsett.counter = 0;
	
	gphsett.resize(2);
	
	gphsett[0].SetGrid(-50, -50, 50, 50, 400, 400);
	gphsett[0].SetPlot(400, 400, 40, 40, 20, 20, "dens");
	gphsett[0].SetAxes("x [kpc]", 10, "y [kpc]", 10);
	gphsett[0].SetColor(1,1,0,20);
	gphsett[0].CalcRatio();
		
	gphsett[1].SetGrid(-50, -800, 50, 800, 400, 400);
	gphsett[1].SetPlot(400, 400, 40, 40, 20, 20, "x_vx");
	gphsett[1].SetAxes("x [kpc]", 10, "vx [km/s]", 100);
	gphsett[1].SetColor(1,1,0,20);
	gphsett[1].CalcRatio();
	
/*	gphsett[2].SetGrid(-50, -400, 50, 400, 400, 400);
	gphsett[2].SetPlot(400, 400, 40, 40, 20, 20, "x_vz");
	gphsett[2].SetAxes("x [kpc]", 10, "vz [km/s]", 100);
	gphsett[2].SetColor(0,0,1.0,20);
	gphsett[2].CalcRatio();
*/
	ss.str(std::string());
	float escape_vel = sqrt(-2*PotPlum(GMP, separace, PPlumSc)-2*PotPlum(GMS, separace, SPlumSc));
	std::cout << escape_vel*VUNIT << "\n" << GMStoP*escape_vel*VUNIT/(1+GMStoP) << "\n" << escape_vel*VUNIT/(1+GMStoP) << "\n";
	
	ss.str(std::string());
	ss << "#plumer1 " << GMP*GMUNIT << " Msun, scale " << PPlumSc << " kpc" << std::endl;
//	ss << "#30000 sticky - SPlumSc, trim 5*PPlumSc, r: 0.015 pc, Cpara 0.5, Cperp 1.0" << std::endl;
	ss << "#plumer2 " << GMS*GMUNIT << " Msun, scale " << SPlumSc << " kpc" << std::endl;
	ss << "#10000 stars - SPlumSc, trim 5*SPlumSc" << std::endl;
	ss << "#10000 sticky - SPlumSc, trim 10*SPlumSc, r: 0.005 pc, Cpara 0.5, Cperp 1.0" << std::endl;
	ss << "#esc vel: " << escape_vel*VUNIT << ", sec " << GMStoP*escape_vel*VUNIT/(1+GMStoP) << ", prim " << escape_vel*VUNIT/(1+GMStoP) << " (km/s)" << std::endl;
	
	escape_vel *= rel_vel;
	obj->AddPotential(T_PLUM,GMP,PPlumSc,0,0,0,0,0,0);
//	GenPlummerStickyTrim2(obj,30000,GMP,PPlumSc,1000/GMUNIT,0.015,5*PPlumSc);
	//GenPlummerStar(obj,10000,GMP,PPlumSc);
	
	obj->AddGroup();
	obj->AddPotential(T_PLUM,GMS,SPlumSc,0,0,0,0,0,0);
	GenPlummerStarTrim2(obj,10000,GMS,SPlumSc,5*SPlumSc);
	GenPlummerStickyTrim2(obj,10000,GMS,SPlumSc,1000/GMUNIT,0.005,10*SPlumSc);
	
	evo->SetStickyCollParams(0.5,1.0,5);
	
	ss << "#rel vel: " << escape_vel*VUNIT << ", sec " << GMStoP*escape_vel*VUNIT/(1+GMStoP) << ", prim " << escape_vel*VUNIT/(1+GMStoP) << " (km/s)" << std::endl;
	ss << "#separation: " << separace << " kpc" << std::endl;
	ss << "#T [Myr]\tN_coll" << std::endl;
	WriteLine(ss.str());
	
	evo->SetObjects(obj);
	evo->InitLeapFrog();
	
/*	sGraphSetting relaxgph;
	relaxgph.SetGrid(-10, -10, 10, 10, 400, 400);
	relaxgph.SetPlot(400, 400, 40, 40, 20, 20, "relax");
	relaxgph.SetAxes("x [kpc]", 2, "y [kpc]", 2);
	relaxgph.SetColor(0,0,1.0,40);
	relaxgph.CalcRatio();
*/	
	int count=0;
	std::cout << "relaxuji..." << std::endl;
	for (int i=0; i<1000; i++) { //relaxace
/*		if (i%2 == 0) {
			relaxgph.data_mriz.ClearData();
			FillDensityGrid(obj, relaxgph, obj->star.size(), obj->star.size()+obj->sticky.size());
			iKresliGrid(relaxgph,1.0,1.0,0,10);
			
			relaxgph.data_mriz.ClearData();
			FillDensityGrid(obj, relaxgph, 0, obj->star.size());
			KresliGrid(relaxgph);
			
			KresliOsy(relaxgph, obj);
			
			UlozPng(relaxgph, count);
			
			count++;
		}
*/		
		evo->Relax(0);
		evo->Relax(1);
	}
	
	obj->SetGroupOffset(0,0,0,0,-GMStoP*rel_x_vel*escape_vel/(1+GMStoP),-GMStoP*sqrt(1-rel_x_vel*rel_x_vel)*escape_vel/(1+GMStoP),0);
	obj->SetGroupOffset(1,-separace,0,0,rel_x_vel*escape_vel/(1+GMStoP),sqrt(1-rel_x_vel*rel_x_vel)*escape_vel/(1+GMStoP),0);
	
	evo->StickyParticlesOn();
	std::cout << "simuluji..." << std::endl;
}

void cPreset::MergerStickyTrimDraw() {
	obj->SetZeroPoint(obj->potential[0].pos[X],obj->potential[0].pos[Y],0);
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], obj->star.size(), obj->star.size()+obj->sticky.size());
	FillPhaseGrid(obj, gphsett[1], obj->star.size(), obj->star.size()+obj->sticky.size(),0,0);
//	FillPhaseGrid(obj, gphsett[2], obj->star.size(), obj->star.size()+obj->sticky.size(),0,2);	

	for (int i=0;i<gphsett.size();i++) {
		iKresliGrid(gphsett[i],1.0,1.0,0,14);
	}	
	
	
/*	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], 0, obj->star.size());
	FillPhaseGrid(obj, gphsett[1], 0, obj->star.size(),0,0);
//	FillPhaseGrid(obj, gphsett[2], 0, obj->star.size(),0,2);
*/			
	
	
	for (int i=0;i<gphsett.size();i++) {
//		KresliGrid(gphsett[i]);
		if (i==0) KresliKolizeXY(gphsett[i], evo);
		KresliOsy(gphsett[i], obj);
//		KresliPaletu(gphsett[i],gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+2, gphsett[i].edge_size[1], gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+gphsett[i].edge_size[2]-2, gphsett[i].canv_size[Y]+gphsett[i].edge_size[1]);
		UlozPng(gphsett[i], gsett.counter);
	}
}

void cPreset::MergerStickyTrimCtrl() {
	if (obj->potential[1].pos[X] >= obj->potential[0].pos[X] && obj->potential[1].active) {
		obj->SetPotentialVel(0,0,0,0);
		obj->SetPotentialVel(1,0,0,0);
		
		obj->TurnOffPotential(1);
	}
	
	ss.str(std::string());
	ss << obj->time << "\t" << evo->total_num_coll << std::endl;
	WriteLine(ss.str());
	
}

void cPreset::MergerStickyTrimFinal() {
	file.close();
	for (int i=0;i<gphsett.size();i++) gphsett[i].plot.close();
}

//-------------------------------------------------------------------------------------------------------

void cPreset::LineInit() {
	file.open("analyza.txt", std::ios::trunc);
	
	gsett.final_time = 1401;
	gsett.draw_every = 2;
	gsett.save = false;
	
	gsett.counter = 0;
	
	gphsett.resize(1);
	
	gphsett[0].SetGrid(-50, -50, 50, 50, 400, 400);
	gphsett[0].SetPlot(400, 400, 40, 40, 20, 20, "dens");
	gphsett[0].SetAxes("x [kpc]", 10, "y [kpc]", 10);
	gphsett[0].SetColor(1.0,1.0,0.0,150);
	gphsett[0].CalcRatio();
		
	ss.str(std::string());
	
	const float delka = 10.0;
	const float pocet = 10000;
	
	ss.str(std::string());
	ss << "#plumer " << GMP*GMUNIT << " Msun, scale " << PPlumSc << " kpc" << std::endl;
	ss << "#" << pocet << " stars" << std::endl;
				
	obj->AddPotential(T_PLUM,GMP,PPlumSc,0,0,0,0,0,0);
	//GenPlummerStar(obj,0,GMP,PPlumSc);
	
	obj->AddGroup();
	GenStarLine(obj, 1, 0, delka/pocet, 0, pocet);
	obj->SetGroupOffset(-20,-delka/2,0,100,0,0);
	
	ss << "#separation: 20 kpc" << std::endl;
//	ss << "#T [Myr]\tN_coll" << std::endl;
	WriteLine(ss.str());
	
	evo->SetObjects(obj);
	evo->InitLeapFrog();
}

void cPreset::LineDraw() {
	obj->SetZeroPoint(obj->potential[0].pos[X],obj->potential[0].pos[Y],0);
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], 0, obj->star.size());
	FillPhaseGrid(obj, gphsett[1], 0, obj->star.size(),0,0);
	FillPhaseGrid(obj, gphsett[2], 0, obj->star.size(),0,2);
	
	for (int i=0;i<gphsett.size();i++) {
		KresliGrid(gphsett[i]);
		KresliOsy(gphsett[i], obj);
		KresliPaletu(gphsett[i],gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+2, gphsett[i].edge_size[1], gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+gphsett[i].edge_size[2]-2, gphsett[i].canv_size[Y]+gphsett[i].edge_size[1]);
		UlozPng(gphsett[i], gsett.counter);
	}

}

void cPreset::LineCtrl() {
	if (obj->potential[1].pos[X] >= obj->potential[0].pos[X] && obj->potential[1].active) {
		obj->SetPotentialVel(0,0,0,0);
		obj->SetPotentialVel(1,0,0,0);
		
		obj->TurnOffPotential(1);
	}
	
//	ss.str(std::string());
//	ss << obj->time << std::endl;
//	WriteLine(ss.str());
}

void cPreset::LineFinal() {
	file.close();
	for (int i=0;i<gphsett.size();i++) gphsett[i].plot.close();
}

//---------------------------------------------------------------------------------------------

void cPreset::HommDiscInit() {
	file.open("analyza.txt", std::ios::trunc);
	
	gsett.final_time = 2001;
	gsett.draw_every = 2;
	gsett.save = false;
	
	gsett.counter = 0;
	
	gphsett.resize(1);
	
	gphsett[0].SetGrid(-5, -5, 5, 5, 200, 200);
	gphsett[0].SetPlot(401, 401, 40, 40, 20, 20, "dens");
	gphsett[0].SetAxes("x [kpc]", 1, "y [kpc]", 1);
	gphsett[0].SetColor(1.0,1.0,0.0,5);
	gphsett[0].CalcRatio();
		
	ss.str(std::string());
	ss << "#plumer hom disc " << GMS*GMUNIT << " Msun, scale " << SPlumSc << " kpc, radius " << 8*SPlumSc << " kpc" <<  std::endl;
	//ss << "#1000000 stars"  << std::endl;
	ss << "#20000 sticky, radius 0.025, koef: para 0.5, perp 1.0, prah 0 km/s, tloustka 100 pc " << std::endl;
	
	obj->AddPotential(T_PLUM,GMS,SPlumSc,0,0,0,0,0,0);
//	GenPlummHommStarDiscHot(obj,1000000,GMS,SPlumSc,8*SPlumSc,0.1);
	GenPlummHommStickyDisc(obj,30000,GMS,SPlumSc,8*SPlumSc,0.25,1,0.01);
	
	evo->SetStickyCollParams(0.5,0.95,5);
	
	ss << "#GMUNIT " << GMUNIT << ", VUNIT" << VUNIT << std::endl;
	
	evo->SetObjects(obj);
	evo->InitLeapFrog();
	evo->StickyParticlesOn();
	
	for (int i=0; i<1000; i++) evo->Relax(0);
}

void cPreset::HommDiscDraw() {
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], 0, obj->star.size());
	FillDensityGrid(obj, gphsett[0], obj->star.size(), obj->star.size()+obj->sticky.size());
		
	for (int i=0;i<gphsett.size();i++) {
		KresliGrid(gphsett[i]);
		KresliOsy(gphsett[i], obj);
//		KresliPaletu(gphsett[i],gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+2, gphsett[i].edge_size[1], gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+gphsett[i].edge_size[2]-2, gphsett[i].canv_size[Y]+gphsett[i].edge_size[1]);
		UlozPng(gphsett[i], gsett.counter);
	}

}

void cPreset::HommDiscCtrl() {
	ss.str(std::string());
	ss << obj->time << "\t" << evo->total_num_coll << std::endl;
	WriteLine(ss.str());
}

void cPreset::HommDiscFinal() {
	file.close();
	for (int i=0;i<gphsett.size();i++) gphsett[i].plot.close();
}

//-----------------------------------------------------------------------------------------------------------------------------

void cPreset::DiscMergerInit() {
	file.open("analyza.txt", std::ios::trunc);
	
	gsett.final_time = 1401;
	gsett.draw_every = 2;
	gsett.save = false;
	gsett.save_every = 1400;
	
	gsett.counter = 0;
	
	gphsett.resize(1);
	
	gphsett[0].SetGrid(-50, -50, 50, 50, 400, 400);
	gphsett[0].SetPlot(400, 400, 40, 40, 20, 20, "dens");
	gphsett[0].SetAxes("x [kpc]", 10, "y [kpc]", 10);
	gphsett[0].SetColor(1.0,1.0,0.0,150);
	gphsett[0].CalcRatio();
/*			
	gphsett[1].SetGrid(-50, -800, 50, 800, 400, 400);
	gphsett[1].SetPlot(400, 400, 40, 40, 20, 20, "x_vx");
	gphsett[1].SetAxes("x [kpc]", 10, "vx [km/s]", 100);
	gphsett[1].SetColor(1.0,1.0,0.0,20);
	gphsett[1].CalcRatio();
	
	gphsett[2].SetGrid(-50, -400, 50, 400, 400, 400);
	gphsett[2].SetPlot(400, 400, 40, 40, 20, 20, "x_vy");
	gphsett[2].SetAxes("x [kpc]", 10, "vy [km/s]", 100);
	gphsett[2].SetColor(1.0,1.0,0.0,20);
	gphsett[2].CalcRatio();
*/	
//	#define Dini 50
	
	ss.str(std::string());
	float escape_vel = sqrt(-2*PotPlum(GMP, Dini, PPlumSc)-2*PotPlum(GMS, Dini, SPlumSc));
	std::cout << escape_vel*VUNIT << "\n" << GMStoP*escape_vel*VUNIT/(1+GMStoP) << "\n" << escape_vel*VUNIT/(1+GMStoP) << "\n";
	
	ss.str(std::string());
	ss << "#plumer1 " << GMP*GMUNIT << " Msun, scale " << PPlumSc << " kpc" << std::endl;
	ss << "#plumer2 " << GMS*GMUNIT << " Msun, scale " << SPlumSc << " kpc" << std::endl;
	ss << "#100000 stars, trim " << 5*SPlumSc << " kpc, tloustka " << 0.3 << " kpc" << std::endl;
	ss << "#20000 sticky, trim " << 10*SPlumSc << " kpc, tloustka " << 0.1 << " kpc" << std::endl;
//	ss << "#dinfric: 0.5 char sekundaru, 2 char primaru" << std::endl;
	ss << "#esc vel: " << escape_vel*VUNIT << ", sec " << GMStoP*escape_vel*VUNIT/(1+GMStoP) << ", prim " << escape_vel*VUNIT/(1+GMStoP) << " (km/s)" << std::endl;
	
	
	const float rel_vel = 0.9;			// > 1 - unbound, = 1 - parabolic, < 1 - bound
	const float rel_x_vel = 1.0;
				
	escape_vel *= rel_vel;
	obj->AddPotential(T_PLUM,GMP,PPlumSc,0,0,0,0,0,0);
	//GenPlummerStar(obj,0,GMP,PPlumSc);
	obj->SetGroupOffset(0,0,0,-GMStoP*rel_x_vel*escape_vel/(1+GMStoP),-GMStoP*sqrt(1-rel_x_vel*rel_x_vel)*escape_vel/(1+GMStoP),0);
	
	obj->AddGroup();
	obj->AddPotential(T_PLUM,GMS,SPlumSc,0,0,0,0,0,0);
	GenPlummHommStarDisc(obj,3000000,GMS,SPlumSc,5*SPlumSc,0.15);
	//GenPlummHommStarDisc(obj,500000,GMS,SPlumSc,8*SPlumSc,0.05);
	//GenPlummHommStickyDisc(obj,20000,GMS,SPlumSc,10*SPlumSc,0.05,1,0.005);
	//obj->SetGroupRot(90,Y);
	obj->SetGroupOffset(-DTiniS-DTiniP,0,0,rel_x_vel*escape_vel/(1+GMStoP),sqrt(1-rel_x_vel*rel_x_vel)*escape_vel/(1+GMStoP),0);
	
	ss << "#rel vel: " << escape_vel*VUNIT << ", sec " << GMStoP*escape_vel*VUNIT/(1+GMStoP) << ", prim " << escape_vel*VUNIT/(1+GMStoP) << " (km/s)" << std::endl;
	ss << "#separation: " << Dini << " kpc" << std::endl;
	ss << "#T [Myr]\tN_coll" << std::endl;
	WriteLine(ss.str());
	
	evo->SetObjects(obj);
	
	evo->SetStickyCollParams(0.5,1.0,0);
	//evo->StickyParticlesOn();
		
	evo->InitLeapFrog();
}

void cPreset::DiscMergerDraw() {
	obj->SetZeroPoint(obj->potential[0].pos[X],obj->potential[0].pos[Y],0);
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	/*
	FillDensityGrid(obj, gphsett[0], 0, obj->star.size()/2);
	FillPhaseGrid(obj, gphsett[1], 0, obj->star.size()/2,0,0);
	FillPhaseGrid(obj, gphsett[2], 0, obj->star.size()/2,0,1);
	*/
	
	FillDensityGrid(obj, gphsett[0], 0, obj->star.size());
//	FillPhaseGrid(obj, gphsett[1], 0, obj->star.size(),0,0);
//	FillPhaseGrid(obj, gphsett[2], 0, obj->star.size(),0,1);
	
	for (int i=0;i<gphsett.size();i++) KresliGrid(gphsett[i]);
//	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	/*
	FillDensityGrid(obj, gphsett[0], obj->star.size()/2, obj->star.size());
	FillPhaseGrid(obj, gphsett[1], obj->star.size()/2, obj->star.size(),0,0);
	FillPhaseGrid(obj, gphsett[2], obj->star.size()/2, obj->star.size(),0,1);
	*/
	
//	FillDensityGrid(obj, gphsett[0], obj->star.size(), obj->sticky.size());
//	FillPhaseGrid(obj, gphsett[1], obj->star.size(), obj->sticky.size(),0,0);
//	FillPhaseGrid(obj, gphsett[2], obj->star.size(), obj->sticky.size(),0,1);
	
	for (int i=0;i<gphsett.size();i++) {
//		iKresliGrid(gphsett[i],0.0,0.0,1.0,40);
		KresliOsy(gphsett[i], obj);
		//KresliPaletu(gphsett[i],gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+2, gphsett[i].edge_size[1], gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+gphsett[i].edge_size[2]-2, gphsett[i].canv_size[Y]+gphsett[i].edge_size[1]);
		UlozPng(gphsett[i], gsett.counter);
	}

}

void cPreset::DiscMergerCtrl() {
	if (obj->potential[1].pos[X] >= obj->potential[0].pos[X] && obj->potential[1].active) {
		obj->SetPotentialVel(0,0,0,0);
		obj->SetPotentialVel(1,0,0,0);
		
		obj->TurnOffPotential(1);
	}
	
/*
	float acc[3];
	
	if (obj->potential[1].active) {
		HandleDynamicalFrictionPlummer(obj->potential[0], obj->potential[1], acc);
		for (int i=0;i<3;i++) {
			if (abs(obj->potential[1].vel[i])<abs(acc[i]) && obj->potential[1].vel[i]*acc[i]<0) obj->potential[1].vel[i]=0;
			else obj->potential[1].vel[i]+=acc[i]*1;
			//std::cout << acc[i] << "	";
		}
	//	std::cout << std::endl;
	}
*/
	
	ss.str(std::string());
	ss << obj->time << "\t" << evo->total_num_coll << std::endl;
	WriteLine(ss.str());
}

void cPreset::DiscMergerFinal() {
	file.close();
	for (int i=0;i<gphsett.size();i++) gphsett[i].plot.close();
}

//--------------------------------------------------------------------------------------------------

void cPreset::DynFricTestInit() {
	file.open("analyza.txt", std::ios::trunc);
	
	gsett.final_time = 5001;
	gsett.draw_every = 2;
	gsett.save = false;
	
	gsett.counter = 0;
	
	gphsett.resize(2);
	
	gphsett[0].SetGrid(-150, -150, 150, 150, 600, 600);
	gphsett[0].SetPlot(600, 600, 40, 40, 20, 20, "dens");
	gphsett[0].SetAxes("x [kpc]", 25, "y [kpc]", 25);
	//gphsett[0].SetColor(1.0,1.0,0.0,100);
	gphsett[0].SetColor(1.0,1.0,0.0,400);
	gphsett[0].CalcRatio();
	
	gphsett[1].SetGrid(-50, -50, 50, 50, 600, 600);
	gphsett[1].SetPlot(600, 600, 40, 40, 20, 20, "dens_zoom");
	gphsett[1].SetAxes("x [kpc]", 10, "y [kpc]", 10);
	//gphsett[0].SetColor(1.0,1.0,0.0,100);
	gphsett[1].SetColor(1.0,1.0,0.0,100);
	gphsett[1].CalcRatio();

	ss.str(std::string());
	ss << "#plumer1 " << GMP*GMUNIT << " Msun, scale " << PPlumSc << " kpc" << std::endl;
	ss << "#plumer2 " << GMS*GMUNIT << " Msun, scale " << SPlumSc << " kpc" << std::endl;
	
	obj->AddPotential(T_PLUM,GMP,PPlumSc,0,0,0,0,0,0);
	//obj->AddPotential(T_PLUM,10*GMP,10*PPlumSc,0,0,0,0,0,0);
	//GenPlummerStar(obj,1000000,GMP,PPlumSc);
	//GenPlummerStar(obj,1000000,10*GMP,10*PPlumSc);
	obj->SetGroupOffset(0,0,0,0,0,0);
	
	obj->AddGroup();
	obj->AddPotential(T_PLUM,GMS,SPlumSc,0,0,0,0,0,0);
	//obj->AddPotential(T_PLUM,10*GMS,10*SPlumSc,0,0,0,0,0,0);
	GenPlummerStar(obj,2000000,GMS,SPlumSc);
	//GenPlummerStar(obj,3000000,10*GMS,10*SPlumSc);
	obj->SetGroupOffset(150,0,0,-0.8*VescPlum(GMS,150,PPlumSc),0,0);
	//obj->SetGroupOffset(0,150,0,sqrt(fabs(150*AccPlum(GMP,150,150*150,PPlumSc*PPlumSc))),0,0);
	//obj->SetGroupOffset(0,150,0,sqrt(fabs(150*AccPlum(10*GMP,150,150*150,10*10*PPlumSc*PPlumSc))),0,0);
	//obj->SetGroupOffset(0,150,0,sqrt(fabs(150*AccPlum(GMP,150,150*150,PPlumSc*PPlumSc)+150*AccPlum(10*GMP,150,150*150,100*PPlumSc*PPlumSc))),0,0);
	//std::cout << sqrt(fabs(10.0*AccPlum(GMP,10,10*10,PPlumSc*PPlumSc)+10*AccPlum(10.0*GMP,10,10*10,100*PPlumSc*PPlumSc))) << std::endl;
	WriteLine(ss.str());
	
	evo->SetObjects(obj);
	evo->InitLeapFrog();
}

void cPreset::DynFricTestDraw() {
	obj->SetZeroPoint(obj->potential[0].pos[X],obj->potential[0].pos[Y],0);
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], 0, obj->star.size());
	FillDensityGrid(obj, gphsett[1], 0, obj->star.size());
	
	for (int i=0;i<gphsett.size();i++) {
		KresliGrid(gphsett[i]);
		/*
		KresliPotencial(gphsett[0], obj->potential[0].pos[X], obj->potential[0].pos[Y], 0.65*sqrt(obj->potential[0].chr2), 0.08, 1.0, 1.0, 0.0);
		KresliPotencial(gphsett[0], obj->potential[0].pos[X], obj->potential[0].pos[Y], sqrt(obj->potential[0].chr2), 0.08, 1.0, 1.0, 0.0);
		KresliPotencial(gphsett[0], obj->potential[0].pos[X], obj->potential[0].pos[Y], 1.3*sqrt(obj->potential[0].chr2), 0.08, 1.0, 1.0, 0.0);
		
		KresliPotencial(gphsett[0], obj->potential[1].pos[X], obj->potential[1].pos[Y], 0.65*sqrt(obj->potential[1].chr2), 0.08, 0.2, 0.2, 1.0);
		KresliPotencial(gphsett[0], obj->potential[1].pos[X], obj->potential[1].pos[Y], sqrt(obj->potential[1].chr2), 0.08, 0.2, 0.2, 1.0);
		KresliPotencial(gphsett[0], obj->potential[1].pos[X], obj->potential[1].pos[Y], 1.3*sqrt(obj->potential[1].chr2), 0.08, 0.2, 0.2, 1.0);
		
		KresliPotencial(gphsett[0], obj->potential[2].pos[X], obj->potential[2].pos[Y], 0.65*sqrt(obj->potential[2].chr2), 0.08, 1.0, 1.0, 0.0);
		KresliPotencial(gphsett[0], obj->potential[2].pos[X], obj->potential[2].pos[Y], sqrt(obj->potential[2].chr2), 0.08, 1.0, 1.0, 0.0);
		KresliPotencial(gphsett[0], obj->potential[2].pos[X], obj->potential[2].pos[Y], 1.3*sqrt(obj->potential[2].chr2), 0.08, 1.0, 1.0, 0.0);
		
		KresliPotencial(gphsett[0], obj->potential[3].pos[X], obj->potential[3].pos[Y], 0.65*sqrt(obj->potential[3].chr2), 0.08, 0.2, 0.2, 1.0);
		KresliPotencial(gphsett[0], obj->potential[3].pos[X], obj->potential[3].pos[Y], sqrt(obj->potential[3].chr2), 0.08, 0.2, 0.2, 1.0);
		KresliPotencial(gphsett[0], obj->potential[3].pos[X], obj->potential[3].pos[Y], 1.3*sqrt(obj->potential[3].chr2), 0.08, 0.2, 0.2, 1.0);
		*/
		KresliOsy(gphsett[i], obj);
//		KresliPaletu(gphsett[i],gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+2, gphsett[i].edge_size[1], gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+gphsett[i].edge_size[2]-2, gphsett[i].canv_size[Y]+gphsett[i].edge_size[1]);
		UlozPng(gphsett[i], gsett.counter);
	}

}

void cPreset::DynFricTestCtrl() {

	float acc[3]; for(int i=0;i<3;i++) acc[i]=0;
	
	HandleDynamicalFrictionPlummer(obj->potential[0], obj->potential[1], acc);
	
	for (int i=0;i<3;i++) obj->potential[1].vel[i]+=acc[i]*DTMyr;
	
	//HandleDynamicalFrictionPlummer(obj->potential[0], obj->potential[2], acc);
	//HandleDynamicalFrictionPlummer(obj->potential[1], obj->potential[2], acc);
	
	//for (int i=0;i<3;i++) obj->potential[2].vel[i]+=acc[i]*DTMyr;
	
	//HandleDynamicalFrictionPlummer(obj->potential[0], obj->potential[3], acc);
	//HandleDynamicalFrictionPlummer(obj->potential[1], obj->potential[3], acc);
	//HandleDynamicalFrictionPlummer(obj->potential[2], obj->potential[3], acc);
	
	//for (int i=0;i<3;i++) obj->potential[3].vel[i]+=acc[i]*DTMyr;
	
	float old_mass, mass_loss;
	old_mass = obj->potential[1].gmass;
	
	mass_loss = OutMassPlum(TidalPlum(obj->potential[0], obj->potential[1]), obj->potential[1])/(TidalPlum(obj->potential[0], obj->potential[1])/DispPlum(obj->potential[1].gmass,0,sqrt(obj->potential[1].chr2)));
	
	obj->potential[1].gmass -= mass_loss*DTMyr;
	obj->potential[1].chr2 = obj->potential[1].chr2 * pow(obj->potential[1].gmass/old_mass,2/3.0);
	
	/*
	float old_mass, mass_loss, mass_loss1, mass_loss2;
	old_mass = obj->potential[2].gmass;
	
	mass_loss1 = OutMassPlum(TidalPlum(obj->potential[0], obj->potential[2]), obj->potential[2])/(TidalPlum(obj->potential[0], obj->potential[2])/DispPlum(obj->potential[2].gmass,0,sqrt(obj->potential[2].chr2)));
	mass_loss2 = OutMassPlum(TidalPlum(obj->potential[1], obj->potential[2]), obj->potential[2])/(TidalPlum(obj->potential[1], obj->potential[2])/DispPlum(obj->potential[2].gmass,0,sqrt(obj->potential[2].chr2)));
	
	if (mass_loss1 > mass_loss2) mass_loss = mass_loss1; else mass_loss = mass_loss2;
	
	obj->potential[2].gmass -= mass_loss*DTMyr;
	obj->potential[2].chr2 = obj->potential[2].chr2 * pow(obj->potential[2].gmass/old_mass,2/3.0);
	*/
	
	/*
	float tidal1, tidal2, mass_loss1, mass_loss2, old_mass;
	
	tidal1 = TidalPlum(obj->potential[0], obj->potential[2]);
	tidal2 = TidalPlum(obj->potential[0], obj->potential[3]);
	
	mass_loss1 = OutMassPlum(tidal1, obj->potential[2])/(tidal1/DispPlum(obj->potential[2].gmass,0,sqrt(obj->potential[2].chr2)));
	mass_loss2 = OutMassPlum(tidal2, obj->potential[3])/(tidal1/DispPlum(obj->potential[3].gmass,0,sqrt(obj->potential[3].chr2)));
	
	tidal1 = TidalPlum(obj->potential[1], obj->potential[2]);
	tidal2 = TidalPlum(obj->potential[1], obj->potential[3]);
	
	mass_loss1 += OutMassPlum(tidal1, obj->potential[2])/(tidal1/DispPlum(obj->potential[2].gmass,0,sqrt(obj->potential[2].chr2)));
	mass_loss2 += OutMassPlum(tidal2, obj->potential[3])/(tidal1/DispPlum(obj->potential[3].gmass,0,sqrt(obj->potential[3].chr2)));
	
	old_mass = obj->potential[2].gmass;
	obj->potential[2].gmass -= mass_loss1*DTMyr;
	obj->potential[2].chr2 = obj->potential[2].chr2 * pow(obj->potential[2].gmass/old_mass,2/3.0);
	
	old_mass = obj->potential[3].gmass;
	obj->potential[3].gmass -= mass_loss2*DTMyr;
	obj->potential[3].chr2 = obj->potential[3].chr2 * pow(obj->potential[3].gmass/old_mass,2/3.0);
	*/
	
	ss.str(std::string());
	ss << obj->time << "	" << fabs(obj->potential[1].pos[X])  << "	" << obj->potential[1].gmass << "	" << sqrt(obj->potential[1].chr2) << std::endl;
	WriteLine(ss.str());
}

void cPreset::DynFricTestFinal() {
	file.close();
	for (int i=0;i<gphsett.size();i++) gphsett[i].plot.close();
}

//----------------------------------------------------------------------------------------------

void cPreset::RozetaInit() {
	file.open("analyza.txt", std::ios::trunc);
	
	gsett.final_time = 1001;
	gsett.draw_every = 2000;
	gsett.save = false;
	gsett.counter = 0;
	
	gphsett.resize(1);
	gphsett[0].SetGrid(-2, -2, 2, 2, 200, 200);
	gphsett[0].SetPlot(201, 201, 40, 40, 20, 20, "dens");
	gphsett[0].SetAxes("x [kpc]", 0.5, "y [kpc]", 0.5);
	gphsett[0].SetColor(1.0,1.0,0.0,1);
	gphsett[0].CalcRatio();
		
	ss.str(std::string());
	ss << "#Plummer " << GMS*GMUNIT << " M_sun, scale " << SPlumSc << " kpc" << std::endl;
	ss << "#3 stars, trimm2 " << 5*SPlumSc  << " plummer, scale " << SPlumSc  << std::endl;
	ss << "#GMUNIT [M_sun] " << GMUNIT << ", VUNIT [km/s]" << VUNIT << std::endl;
	
	ss << "#T [Myr]";
	for (int i=1; i<=3; i++) ss << "\t x" << i << "\t y" << i;
	ss << std::endl;
	WriteLine(ss.str());
	
	obj->AddPotential(T_PLUM,GMS,SPlumSc,0,0,0,0,0,0);
	obj->AddStar(1,-2,0,0,0,0.05,0);
	obj->AddStar(1,-1,0,0,0,0.05,0);
	obj->AddStar(1,-0.5,0,0,0,0.05,0);
	
	evo->SetObjects(obj);
	evo->InitLeapFrog();
	evo->StickyParticlesOff();
}

void cPreset::RozetaDraw() {
	
	for (int i=0;i<gphsett.size();i++) gphsett[i].data_mriz.ClearData();
	
	FillDensityGrid(obj, gphsett[0], 0, obj->star.size());
		
	for (int i=0;i<gphsett.size();i++) {
		KresliGrid(gphsett[i]);
		KresliOsy(gphsett[i], obj);
//		KresliPaletu(gphsett[i],gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+2, gphsett[i].edge_size[1], gphsett[i].canv_size[X]+gphsett[i].edge_size[0]+gphsett[i].edge_size[2]-2, gphsett[i].canv_size[Y]+gphsett[i].edge_size[1]);
		UlozPng(gphsett[i], gsett.counter);
	}

}

void cPreset::RozetaCtrl() {

	ss.str(std::string());
	ss << obj->time ;
	for (int i=0; i<3; i++) ss <<"\t" << obj->star[i].pos[X] << "\t" << obj->star[i].pos[Y];
	ss << std::endl;
	
	WriteLine(ss.str());
}

void cPreset::RozetaFinal() {
	file.close();
	for (int i=0;i<gphsett.size();i++) gphsett[i].plot.close();
}
