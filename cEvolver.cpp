#include "cEvolver.h"

cEvolver::cEvolver() {
	test_particles = true;
	
	Reset();
};

void cEvolver::Reset() {
	dt = DTMyr;
	
	acc_obj.clear();
	acc_pot.clear();
	objects = NULL;
	total_num_coll = 0;
	
	acc_thread.resize(NUM_THREADS);
	coll_thread.resize(NUM_THREADS);
	num_coll_thread.resize(NUM_THREADS);
};

void cEvolver::ResetAcc() {
	#pragma omp for
	for (int i=0;i<acc_obj.size();i++) for (int j=0;j<DIM;j++) acc_obj[i].acc[j] = 0;
	
	#pragma omp for
	for (int i=0;i<acc_pot.size();i++) for (int j=0;j<DIM;j++) acc_pot[i].acc[j] = 0;
};

void cEvolver::SetTimeStep(float a) {
	if (a > 0) dt = a;
};

void cEvolver::SetObjects(cObjects* a) {
	objects = a;
	acc_obj.resize(objects->sg_size);
	acc_pot.resize(objects->potential.size());
};

void cEvolver::UpdateAcc() {
	float R2, temp;
	int thr_n = omp_get_thread_num();
	
	
	std::vector<S_Acc > localAcc(acc_obj.size());
	acc_thread[thr_n] = localAcc.data();
	
	
	if (objects != NULL) {
		
		ResetAcc();
		
		#pragma omp for
		for (int i=0;i<acc_obj.size();i++) {
						
			// výpočet zrychlení interakce potenciál-objekt
			for (int j=0;j<acc_pot.size();j++) if (objects->potential[j].active) {
				R2 = Rad2(objects->sg_object(i)->pos,objects->potential[j].pos);
				for (int k=X;k<DIM;k++) acc_obj[i].acc[k] += AccPotential(objects->potential[j].type, objects->potential[j].gmass, objects->sg_object(i)->pos[k]-objects->potential[j].pos[k], R2, objects->potential[j].chr2);
			}
		
			if (!test_particles) for (int j=i+1;j<acc_obj.size();j++) {
				R2 = Rad2(objects->sg_object(i)->pos,objects->sg_object(j)->pos);
				for (int k=X;k<DIM;k++) {
					
					temp = AccNewt(1, objects->sg_object(i)->pos[k]-objects->sg_object(j)->pos[k], R2);
			
					localAcc[i].acc[k] += temp*objects->sg_object(j)->gmass;
					localAcc[j].acc[k] -= temp*objects->sg_object(i)->gmass;
				}
			}
							
		}
		
		#pragma omp for
		for (int i=0;i<acc_obj.size();i++) for (int j=0;j<NUM_THREADS;j++) for (int k=X;k<DIM;k++) acc_obj[i].acc[k] += acc_thread[j][i].acc[k];
		
		#pragma omp single
		{
		    // výpočet zrychlení interakce potenciál-potenciál
		    for (int i=0;i<acc_pot.size();i++) if (objects->potential[i].active) {
		    	for (int j=i+1;j<acc_pot.size();j++) if (objects->potential[j].active) {
		    		R2 = Rad2(objects->potential[i].pos,objects->potential[j].pos);
		    		for (int k=X;k<DIM;k++) {		 
					acc_pot[i].acc[k] += AccPotential(objects->potential[j].type, objects->potential[j].gmass, objects->potential[i].pos[k]-objects->potential[j].pos[k], R2, objects->potential[j].chr2);
					acc_pot[j].acc[k] += AccPotential(objects->potential[i].type, objects->potential[i].gmass, objects->potential[j].pos[k]-objects->potential[i].pos[k], R2, objects->potential[i].chr2);
				}
			}
		    }
		}
	}
};

void cEvolver::UpdateVel() {	
	    #pragma omp for
	    for (int i=0;i<acc_obj.size();i++) for (int j=X;j<DIM;j++) objects->sg_object(i)->vel[j] += acc_obj[i].acc[j] * dt;
	    
	    #pragma omp for
	    for (int i=0;i<acc_pot.size();i++) if (objects->potential[i].active) for (int j=X;j<DIM;j++) objects->potential[i].vel[j] += acc_pot[i].acc[j] * dt;
};

void cEvolver::UpdatePos() {	
	    #pragma omp for
	    for (int i=0;i<acc_obj.size();i++) for (int j=0;j<DIM;j++) objects->sg_object(i)->pos[j] += objects->sg_object(i)->vel[j] * dt;
	
	    #pragma omp for
	    for (int i=0;i<acc_pot.size();i++) if (objects->potential[i].active) for (int j=0;j<DIM;j++) objects->potential[i].pos[j] += objects->potential[i].vel[j] * dt;
};

void cEvolver::SetStickyCollParams(float para, float perp, float min_vel) {
	coll_para = para;
	coll_perp = perp;
	coll_min_vel2 = min_vel*min_vel;
};

void cEvolver::HandleStickyCollision() {
	float dist, rel_vel2;
	int thr_n = omp_get_thread_num();
	
	S_StColl temp;
	std::vector<S_StColl > localColl;
	
	
		
	#pragma omp for
	for (int i=0;i<objects->sticky.size();i++) for (int j=i+1;j<objects->sticky.size();j++) if (fabs(objects->sticky[i].pos[X]-objects->sticky[j].pos[X])<=(objects->sticky[i].R + objects->sticky[j].R) || fabs(objects->sticky[i].pos[Y]-objects->sticky[j].pos[Y])<=(objects->sticky[i].R + objects->sticky[j].R) || fabs(objects->sticky[i].pos[Z]-objects->sticky[j].pos[Z])<=(objects->sticky[i].R + objects->sticky[j].R)) {
		rel_vel2 = Rad2(objects->sticky[i].vel,objects->sticky[j].vel)*VUNIT*VUNIT;
		if (rel_vel2 > coll_min_vel2) {
			dist = sqrt( Rad2(objects->sticky[i].pos,objects->sticky[j].pos) );
			
			if (dist <= objects->sticky[i].R + objects->sticky[j].R) {
				temp.st1 = &objects->sticky[i];
				temp.st2 = &objects->sticky[j];
				temp.dist = dist;
			
				localColl.push_back(temp);
			}
		}
	}
	
	num_coll_thread[thr_n] = localColl.size();
	coll_thread[thr_n] = localColl.data();		//MUSÍ SE PŘIŘAZOVAT ZDE, NE PŘED FOR, PŘI push_back SE MENI UKAZATEL DAT !!!
	
	#pragma omp barrier
	
	#pragma omp single
	{
	S_Pos kde;
	
	total_num_coll = 0;
	for (int i=0;i<NUM_THREADS;i++)	total_num_coll += num_coll_thread[i];
	
	coll_pos.clear();
	for (int i=0;i<NUM_THREADS;i++) for (int j=0;j<num_coll_thread[i];j++) {
		HandleCollision(coll_thread[i][j], coll_para, coll_perp);
		
		for (int k=0; k<3; k++) kde.pos[k] = (coll_thread[i][j].st1->pos[k]+coll_thread[i][j].st2->pos[k])/2;
		coll_pos.push_back(kde);
	}
	}
};

void cEvolver::InitLeapFrog() {
	#pragma omp parallel num_threads(NUM_THREADS)
	{
	    UpdateAcc();
	    
	    #pragma omp for
	    for (int i=0;i<acc_obj.size();i++) for (int j=X;j<DIM;j++) objects->sg_object(i)->vel[j] -= acc_obj[i].acc[j] * dt/2;
	    
	    #pragma omp for
	    for (int i=0;i<acc_pot.size();i++) if (objects->potential[i].active) for (int j=0;j<DIM;j++) objects->potential[i].vel[j] -= acc_pot[i].acc[j] * dt/2;
	}
};

void cEvolver::Evolve() {
	#pragma omp parallel num_threads(NUM_THREADS)
	{
	    UpdateAcc();
	    UpdateVel();		
	    UpdatePos();

	    if (sticky_particles) HandleStickyCollision();
	}
	objects->time += dt;
};

//------------------------------------------------

void cEvolver::UpdateAcc(int G) {
	float R2, temp;
	int thr_n = omp_get_thread_num();
	
	
	std::vector<S_Acc > localAcc(objects->group[G].NStot);
	acc_thread[thr_n] = localAcc.data();
	
	
	if (objects != NULL) {
		
		ResetAcc();
		
		#pragma omp for
		for (int i=0;i<objects->group[G].NStot;i++) {
						
			// výpočet zrychlení interakce potenciál-objekt
			for (int j=objects->group[G].LPots;j<objects->group[G].LPots+objects->group[G].NPots;j++) if (objects->potential[j].active) {
				R2 = Rad2(objects->sg_object(G,i)->pos,objects->potential[j].pos);
				for (int k=X;k<DIM;k++) acc_obj[i].acc[k] += AccPotential(objects->potential[j].type, objects->potential[j].gmass, objects->sg_object(G,i)->pos[k]-objects->potential[j].pos[k], R2, objects->potential[j].chr2);
			}
		
			if (!test_particles) for (int j=i+1;j<objects->group[G].NStot;j++) {
				R2 = Rad2(objects->sg_object(G,i)->pos,objects->sg_object(G,j)->pos);
				for (int k=X;k<DIM;k++) {
					
					temp = AccNewt(1, objects->sg_object(G,i)->pos[k]-objects->sg_object(G,j)->pos[k], R2);
			
					localAcc[i].acc[k] += temp*objects->sg_object(G,j)->gmass;
					localAcc[j].acc[k] -= temp*objects->sg_object(G,i)->gmass;
				}
			}
							
		}
		
		#pragma omp for
		for (int i=0;i<objects->group[G].NStot;i++) for (int j=0;j<NUM_THREADS;j++) for (int k=X;k<DIM;k++) acc_obj[i].acc[k] += acc_thread[j][i].acc[k];
		
		#pragma omp single
		{
		    // výpočet zrychlení interakce potenciál-potenciál
		    for (int i=objects->group[G].LPots;i<objects->group[G].LPots+objects->group[G].NPots;i++) if (objects->potential[i].active) {
		    	for (int j=i+1;j<objects->group[G].LPots+objects->group[G].NPots;j++) if (objects->potential[j].active) {
		    		R2 = Rad2(objects->potential[i].pos,objects->potential[j].pos);
		    		for (int k=X;k<DIM;k++) {		 
					acc_pot[i].acc[k] += AccPotential(objects->potential[j].type, objects->potential[j].gmass, objects->potential[i].pos[k]-objects->potential[j].pos[k], R2, objects->potential[j].chr2);
					acc_pot[j].acc[k] += AccPotential(objects->potential[i].type, objects->potential[i].gmass, objects->potential[j].pos[k]-objects->potential[i].pos[k], R2, objects->potential[i].chr2);
				}
			}
		    }
		}
	}
};

void cEvolver::UpdateVel(int G) {	
	    #pragma omp for
	    for (int i=0;i<objects->group[G].NStot;i++) for (int j=X;j<DIM;j++) objects->sg_object(G,i)->vel[j] += acc_obj[i].acc[j] * dt;
	    
	    #pragma omp for
	    for (int i=objects->group[G].LPots;i<objects->group[G].LPots+objects->group[G].NPots;i++) if (objects->potential[i].active) for (int j=X;j<DIM;j++) objects->potential[i].vel[j] += acc_pot[i].acc[j] * dt;
};

void cEvolver::UpdatePos(int G) {	
	    #pragma omp for
	    for (int i=0;i<objects->group[G].NStot;i++) for (int j=0;j<DIM;j++) objects->sg_object(G,i)->pos[j] += objects->sg_object(G,i)->vel[j] * dt;
	
	    #pragma omp for
	    for (int i=objects->group[G].LPots;i<objects->group[G].LPots+objects->group[G].NPots;i++) if (objects->potential[i].active) for (int j=0;j<DIM;j++) objects->potential[i].pos[j] += objects->potential[i].vel[j] * dt;
};

void cEvolver::Relax(int G) {
	#pragma omp parallel num_threads(NUM_THREADS)
	{
	    UpdateAcc(G);
	    UpdateVel(G);	
	    UpdatePos(G);
	}
	objects->group[G].relax_time += dt;
};