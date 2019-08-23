#include "all.h"

using namespace std;

// Ray class functions:
int Ray::clear_ray (){
	
	delete []dt_u;
	delete []Iu;
	//cout << "Ray deallocated" << endl;
	return 0;
}

// Ray manipulation:

int Point::allocate_rays (){
	
	int NT = Globals::getNT();
	int g, m;
	
	rays = new (nothrow) Ray * [NT];
	for (g=0; g<NT; g++)
		rays[g] = new (nothrow) Ray [Globals::getNP_C(g)];
		
	return 0;
}

int Point::deallocate_rays (){
	
	for (int g=0; g<Globals::getNT(); g++)
		delete []rays[g];
	delete []rays;
}

// Computation of the intersection points and allocation/deallocation/computation of transfer parameters

int Point::compute_upwind_points (int g, int m, Point ** Grid){
	
	// Needed helping variables:
	
	double deltaX;
	double deltaY;
	
	double deltaXmax;
	double deltaYmax;
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	// There are in total four options where the ray is oriented, and subsequently where is the upwind point:
	
	// Option ONE - LEFT, DOWN
	
	if (Globals::getphi_C(g, m) <= pi/2){
		
		if (i == 0 || j == 0){ // If it is on the edge - no upwind
			rays[g][m].Xu = X;
			rays[g][m].Yu = Y;
			rays[g][m].type_u = 1;
		}
		
		else {
			deltaX = X - Grid[i-1][j].X;
			deltaY = Y - Grid[i][j-1].Y;
			
			deltaYmax = deltaX * Globals::getsinphi_C(g, m) / Globals::getcosphi_C(g,m);
			
			if (deltaYmax <= deltaY){
				deltaXmax = deltaX;
				rays[g][m].type_u = 2;
			}
				
			else {
				deltaYmax = deltaY;
				deltaXmax = deltaY * Globals::getcosphi_C(g, m) / Globals::getsinphi_C(g, m);
				rays[g][m].type_u = 1;
			}
			
			rays[g][m].Xu = X - deltaXmax;
			rays[g][m].Yu = Y - deltaYmax;
		}

	}
	
	// Option TWO - RIGHT, DOWN
	
	else if (Globals::getphi_C(g, m) <= pi){
		
		if (i == NX-1 || j == 0){ // If it is on the edge - no upwind
			rays[g][m].Xu = X;
			rays[g][m].Yu = Y;
			rays[g][m].type_u = 1;
		}
		
		else {
			deltaX = Grid[i+1][j].X - X;
			deltaY = Y - Grid[i][j-1].Y;
			
			deltaYmax = - deltaX * Globals::getsinphi_C(g, m) / Globals::getcosphi_C(g,m);
			
			if (deltaYmax <= deltaY){
				deltaXmax = deltaX;
				rays[g][m].type_u = 2;
			}
				
			else {
				deltaYmax = deltaY;
				deltaXmax = - deltaY * Globals::getcosphi_C(g, m) / Globals::getsinphi_C(g, m);
				rays[g][m].type_u = 1;
			}
			
			rays[g][m].Xu = X + deltaXmax;
			rays[g][m].Yu = Y - deltaYmax;
		}
	}
	
	// Option THREE - RIGHT, UP
	
	else if (Globals::getphi_C(g, m) <= 3 * pi / 2){
		
		if (i == NX-1 || j == NY-1){ // If it is on the edge - no upwind
			rays[g][m].Xu = X;
			rays[g][m].Yu = Y;
			rays[g][m].type_u = 1;
		}
		
		else {
			deltaX = Grid[i+1][j].X - X;
			deltaY = Grid[i][j+1].Y - Y;
			
			deltaYmax = deltaX * Globals::getsinphi_C(g, m) / Globals::getcosphi_C(g,m);
			
			if (deltaYmax <= deltaY){
				deltaXmax = deltaX;
				rays[g][m].type_u = 2;
			}
				
			else {
				deltaYmax = deltaY;
				deltaXmax = deltaY * Globals::getcosphi_C(g, m) / Globals::getsinphi_C(g, m);
				rays[g][m].type_u = 1;
			}
			
			rays[g][m].Xu = X + deltaXmax;
			rays[g][m].Yu = Y + deltaYmax;
		}
	}
	
	// Option FOUR- LEFT, UP
	
	else if (Globals::getphi_C(g, m) <= 2 * pi){
		
		if (i == 0 || j == NY-1){ // If it is on the edge - no upwind
			rays[g][m].Xu = X;
			rays[g][m].Yu = Y;
			rays[g][m].type_u = 1;
		}
		
		else {
			deltaX = X - Grid[i-1][j].X;
			deltaY = Grid[i][j+1].Y - Y;
			
			deltaYmax = - deltaX * Globals::getsinphi_C(g, m) / Globals::getcosphi_C(g,m);
			
			if (deltaYmax <= deltaY){
				deltaXmax = deltaX;
				rays[g][m].type_u = 2;
			}
				
			else {
				deltaYmax = deltaY;
				deltaXmax = - deltaY * Globals::getcosphi_C(g, m) / Globals::getsinphi_C(g, m);
				rays[g][m].type_u = 1;
			}
			
			rays[g][m].Xu = X - deltaXmax;
			rays[g][m].Yu = Y + deltaYmax;
		}
	}
	
	else {
		cout << "Bad angle Phi" << endl;
		return 1;
	}

	// Modifications for the case of the periodic boundary conditions implemented by the approach 2:

	if (i ==0 && Globals::getcosphi_C(g,m) >= 0 && periodic_boundary == 1){ // This is the ray incoming from the left 
		// Toward the left boundary, it is identical to one at i == NX-1

		Grid[NX-1][j].allocate_rays();
		Grid[NX-1][j].compute_upwind_points(g,m, Grid);

		Grid[i][j].rays[g][m].Xu = Grid[NX-1][j].rays[g][m].Xu - Globals::getX(NX-1);
		Grid[i][j].rays[g][m].Yu = Grid[NX-1][j].rays[g][m].Yu;
		Grid[i][j].rays[g][m].type_u = Grid[NX-1][j].rays[g][m].type_u;

		Grid[NX-1][j].deallocate_rays();  

	}

	if (i ==NX-1 && Globals::getcosphi_C(g,m) <= 0 && periodic_boundary == 1){ // This is the ray incoming from the right 
		// Toward the right boundary, it is identical to one at i == 0

		Grid[0][j].allocate_rays();
		Grid[0][j].compute_upwind_points(g,m, Grid);

		Grid[i][j].rays[g][m].Xu = Grid[0][j].rays[g][m].Xu + Globals::getX(NX-1);
		Grid[i][j].rays[g][m].Yu = Grid[0][j].rays[g][m].Yu;
		Grid[i][j].rays[g][m].type_u = Grid[0][j].rays[g][m].type_u;

		Grid[0][j].deallocate_rays();  

	}
	
	return 0;
}

int Point::allocate_pqr() {
	
	// Here we allocate the memory for transport coefficients:
	
	p = new (nothrow) double ** [Globals::getNT()];
	q = new (nothrow) double ** [Globals::getNT()];
	r = new (nothrow) double ** [Globals::getNT()];
	exxp = new (nothrow) double ** [Globals::getNT()];
	
	for (int g = 0; g<Globals::getNT(); g++){
		
		p[g] = new (nothrow) double * [Globals::getNP_C(g)];
		q[g] = new (nothrow) double * [Globals::getNP_C(g)];
		r[g] = new (nothrow) double * [Globals::getNP_C(g)];
		exxp[g] = new (nothrow) double * [Globals::getNP_C(g)];
		
		for (int m = 0; m<Globals::getNP_C(g); m++){
			
			p[g][m] = new (nothrow) double [Globals::getN()];
			q[g][m] = new (nothrow) double [Globals::getN()];
			r[g][m] = new (nothrow) double [Globals::getN()];
			exxp[g][m] = new (nothrow) double [Globals::getN()];
		}
	}
	
	return 0;
}

int Point::deallocate_pqr() {
	
	// Deallocate memory for transport coefficients:
	
	for (int g=0; g<Globals::getNT(); g++){
		
		for (int m = 0; m<Globals::getNP_C(g); m++){
			
			delete []p[g][m];
			delete []q[g][m];
			delete []r[g][m];
			delete []exxp[g][m];
		}
		
		delete []p[g];
		delete []q[g];
		delete []r[g];
		delete []exxp[g];
	}
	
	delete []p;
	delete []q;
	delete []r;
	delete []exxp;
	
	return 0;
	
}

int Point::allocate_rxry (){
	
	rx = new (nothrow) double ** [Globals::getNT()];
	ry = new (nothrow) double ** [Globals::getNT()];
	
	for (int g = 0; g < Globals::getNT(); g++){
		
		rx[g] = new (nothrow) double * [Globals::getNP_C(g)];
		ry[g] = new (nothrow) double * [Globals::getNP_C(g)];
		
		for (int m = 0; m<Globals::getNP_C(g); m++){
			rx[g][m] = new (nothrow) double [Globals::getN()];
			ry[g][m] = new (nothrow) double [Globals::getN()];
		}
	}
	
	return 0;
}

int Point::deallocate_rxry (){
	
	for (int g = 0; g<Globals::getNT(); g++){
		
		for (int m = 0; m<Globals::getNP_C(g); m++){
			
			delete []rx[g][m];
			delete []ry[g][m];
		}
		
		delete []rx[g];
		delete []ry[g];
	}
	
	delete []rx;
	delete []ry;
		
	
	return 0;
}

int Point::allocate_rs (){
	
	// Here we allocate r1- r4 coefficients:
	
	r1 = new (nothrow) double ** [Globals::getNT()];
	r2 = new (nothrow) double ** [Globals::getNT()];
	r3 = new (nothrow) double ** [Globals::getNT()];
	r4 = new (nothrow) double ** [Globals::getNT()];
	
	for (int g = 0; g<Globals::getNT(); g++){
		
		r1[g] = new (nothrow) double * [Globals::getNP_C(g)];
		r2[g] = new (nothrow) double * [Globals::getNP_C(g)];
		r3[g] = new (nothrow) double * [Globals::getNP_C(g)];
		r4[g] = new (nothrow) double * [Globals::getNP_C(g)];
		
		for (int m = 0; m<Globals::getNP_C(g); m++){
			
			r1[g][m] = new (nothrow) double [Globals::getN()];
			r2[g][m] = new (nothrow) double [Globals::getN()];
			r3[g][m] = new (nothrow) double [Globals::getN()];
			r4[g][m] = new (nothrow) double [Globals::getN()];
		}
	}
	
	return 0;
}

int Point::deallocate_rs () {
	
	for (int g = 0; g<Globals::getNT(); g++){
		
		for (int m =0; m<Globals::getNP_C(g); m++){
			
			delete []r1[g][m];
			delete []r2[g][m];
			delete []r3[g][m];
			delete []r4[g][m];
		}
		
		delete []r1[g];
		delete []r2[g];
		delete []r3[g];
		delete []r4[g];
		
	}
	
	delete []r1;
	delete []r2;
	delete []r3;
	delete []r4;
	
	return 0;
}

int Point::allocate_r_full (){
	
	// Here we allocate r1- r4 coefficients:
	
	r1 = new (nothrow) double ** [Globals::getNT()];
	r2 = new (nothrow) double ** [Globals::getNT()];
	r3 = new (nothrow) double ** [Globals::getNT()];
	r4 = new (nothrow) double ** [Globals::getNT()];
	r5 = new (nothrow) double ** [Globals::getNT()];
	r6 = new (nothrow) double ** [Globals::getNT()];
	r7 = new (nothrow) double ** [Globals::getNT()];
	r8 = new (nothrow) double ** [Globals::getNT()];
	
	for (int g = 0; g<Globals::getNT(); g++){
		
		r1[g] = new (nothrow) double * [Globals::getNP_C(g)];
		r2[g]= new (nothrow) double * [Globals::getNP_C(g)];
		r3[g] = new (nothrow) double * [Globals::getNP_C(g)];
		r4[g] = new (nothrow) double * [Globals::getNP_C(g)];
		r5[g] = new (nothrow) double * [Globals::getNP_C(g)];
		r6[g] = new (nothrow) double * [Globals::getNP_C(g)];
		r7[g] = new (nothrow) double * [Globals::getNP_C(g)];
		r8[g] = new (nothrow) double * [Globals::getNP_C(g)];
		
		for (int m = 0; m<Globals::getNP_C(g); m++){
			
			r1[g][m] = new (nothrow) double [Globals::getN()];
			r2[g][m] = new (nothrow) double [Globals::getN()];
			r3[g][m] = new (nothrow) double [Globals::getN()];
			r4[g][m] = new (nothrow) double [Globals::getN()];
			r5[g][m] = new (nothrow) double [Globals::getN()];
			r6[g][m] = new (nothrow) double [Globals::getN()];
			r7[g][m] = new (nothrow) double [Globals::getN()];
			r8[g][m] = new (nothrow) double [Globals::getN()];
			
			for (int n=0; n<Globals::getN(); n++){
				r1[g][m][n] = 0.0;
				r2[g][m][n] = 0.0;
				r3[g][m][n] = 0.0;
				r4[g][m][n] = 0.0;
				r5[g][m][n] = 0.0;
				r6[g][m][n] = 0.0;
				r7[g][m][n] = 0.0;
				r8[g][m][n] = 0.0;
			}
		}
	}
	
	return 0;
}

int Point::deallocate_r_full (){
	
	for (int g = 0; g<Globals::getNT(); g++){
		
		for (int m =0; m<Globals::getNP_C(g); m++){
			
			delete []r1[g][m];
			delete []r2[g][m];
			delete []r3[g][m];
			delete []r4[g][m];
			delete []r5[g][m];
			delete []r6[g][m];
			delete []r7[g][m];
			delete []r8[g][m];
		}
		
		delete []r1[g];
		delete []r2[g];
		delete []r3[g];
		delete []r4[g];
		delete []r5[g];
		delete []r6[g];
		delete []r7[g];
		delete []r8[g];
		
	}
	
	delete []r1;
	delete []r2;
	delete []r3;
	delete []r4;
	delete []r5;
	delete []r6;
	delete []r7;
	delete []r8;
	
	return 0;
}

// Interpolation:
int Point::interpolate_simple (int g, int m, Point ** Grid){
	
	// This is the function where we interpolate Source function, and, if needed, other scalar quantities:
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	double ksi1, ksi2, ksi3; // Coorinates needed for the interpolation.
	double ksi; // Coordinate to perform an interpolation
	int i1, i2, i3, j1, j2, j3;
	
	int dointerpol = 0;
	
	if (Globals::getphi_C(g, m) <= pi/2){ // FIRST QUADRANT:
	
		if (i == 0 || j == 0){
			
			rays[g][m].Su = 0.0;
			dointerpol = 0;
		}
		
		else if (rays[g][m].type_u ==1) { // Interpolation along X
		
			// This is definite:
			
			ksi = rays[g][m].Xu;
			dointerpol = 1;
		
			if (i == 1){
				i1 = i2 = 0;
				i3 = 1;
				j1 = j2 = j3 = j-1;
				ksi1 = Globals::getX(i1);
				ksi2 = Globals::getX(i2);
				ksi3 = Globals::getX(i3);
			}
			else {
				i1 = i-2;
				i2 = i-1;
				i3 = i;
				j1 = j2 = j3 = j-1;
				ksi1 = Globals::getX(i1);
				ksi2 = Globals::getX(i2);
				ksi3 = Globals::getX(i3);
			}
		
		}
		
		else if (rays[g][m].type_u == 2) { // Interpolation along Y
		
			ksi = rays[g][m].Yu;
			dointerpol = 1;
			
			if (j == 1){
				i1 = i2 = i3 = i-1;
				j1 = j2 = j-1;
				j3 = j;
				ksi1 = Globals::getY(j1);
				ksi2 = Globals::getY(j2);
				ksi3 = Globals::getY(j3);
			} 
			
			else {
				i1 = i2 = i3 = i-1;
				j1 = j-2;
				j2 = j-1;
				j3 = j;
				ksi1 = Globals::getY(j1);
				ksi2 = Globals::getY(j2);
				ksi3 = Globals::getY(j3);
			}
		} 
	}
			
	else if (Globals::getphi_C(g, m) <= pi){ // SECOND QUADRANT:
		if (i == NX-1 || j == 0){
			
			rays[g][m].Su = 0.0;
			dointerpol = 0;
		}
		
		else if (rays[g][m].type_u ==1) { // Interpolation along X
		
			// This is definite:
			
			ksi = rays[g][m].Xu;
			dointerpol = 1;
		
			if (i == NX-2){
				i1 = i2 = NX-2;
				i3 = NX-1;
				j1 = j2 = j3 = j-1;
				ksi1 = Globals::getX(i1);
				ksi2 = Globals::getX(i2);
				ksi3 = Globals::getX(i3);
			}
			else {
				i1 = i;
				i2 = i+1;
				i3 = i+2;
				j1 = j2 = j3 = j-1;
				ksi1 = Globals::getX(i1);
				ksi2 = Globals::getX(i2);
				ksi3 = Globals::getX(i3);
			}
		
		}
		
		else if (rays[g][m].type_u == 2) { // Interpolation along Y
		
			ksi = rays[g][m].Yu;
			dointerpol = 1;
			
			if (j == 1){
				i1 = i2 = i3 = i+1;
				j1 = j2 = j-1;
				j3 = j;
				ksi1 = Globals::getY(j1);
				ksi2 = Globals::getY(j2);
				ksi3 = Globals::getY(j3);
			} 
			
			else {
				i1 = i2 = i3 = i+1;
				j1 = j-2;
				j2 = j-1;
				j3 = j;
				ksi1 = Globals::getY(j1);
				ksi2 = Globals::getY(j2);
				ksi3 = Globals::getY(j3);
			}
		}
	}
	
	else if (Globals::getphi_C(g, m) <= pi * 1.5){ // THIRD QUADRANT:
		if (i == NX-1 || j == NY-1){
			
			rays[g][m].Su = 0.0;
			dointerpol = 0;
		}
		
		else if (rays[g][m].type_u ==1) { // Interpolation along X
		
			// This is definite:
			
			ksi = rays[g][m].Xu;
			dointerpol = 1;
		
			if (i == NX-2){
				i1 = i2 = NX-2;
				i3 = NX-1;
				j1 = j2 = j3 = j+1;
				ksi1 = Globals::getX(i1);
				ksi2 = Globals::getX(i2);
				ksi3 = Globals::getX(i3);
			}
			else {
				i1 = i;
				i2 = i+1;
				i3 = i+2;
				j1 = j2 = j3 = j+1;
				ksi1 = Globals::getX(i1);
				ksi2 = Globals::getX(i2);
				ksi3 = Globals::getX(i3);
			}
		
		}
		
		else if (rays[g][m].type_u == 2) { // Interpolation along Y
		
			ksi = rays[g][m].Yu;
			dointerpol = 1;
			
			if (j == NY-2){
				i1 = i2 = i3 = i+1;
				j1 = j2 = j;
				j3 = j+1;
				ksi1 = Globals::getY(j1);
				ksi2 = Globals::getY(j2);
				ksi3 = Globals::getY(j3);
			} 
			
			else {
				i1 = i2 = i3 = i+1;
				j1 = j;
				j2 = j+1;
				j3 = j+2;
				ksi1 = Globals::getY(j1);
				ksi2 = Globals::getY(j2);
				ksi3 = Globals::getY(j3);
			}
		}
	}
	
	else if (Globals::getphi_C(g, m) <= pi * 2){ // FOURTH QUADRANT:
		if (i == 0 || j == NY-1){
			
			rays[g][m].Su = 0.0;
			dointerpol = 0;
		}
		
		else if (rays[g][m].type_u ==1) { // Interpolation along X
		
			// This is definite:
			
			ksi = rays[g][m].Xu;
			dointerpol = 1;
		
			if (i == 1){
				i1 = i2 = 0;
				i3 = 1;
				j1 = j2 = j3 = j+1;
				ksi1 = Globals::getX(i1);
				ksi2 = Globals::getX(i2);
				ksi3 = Globals::getX(i3);
			}
			else {
				i1 = i-2;
				i2 = i-1;
				i3 = i;
				j1 = j2 = j3 = j+1;
				ksi1 = Globals::getX(i1);
				ksi2 = Globals::getX(i2);
				ksi3 = Globals::getX(i3);
			}
		
		}
		
		else if (rays[g][m].type_u == 2) { // Interpolation along Y
		
			ksi = rays[g][m].Yu;
			dointerpol = 1;
			
			if (j == NY-2){
				i1 = i2 = i3 = i-1;
				j1 = j2 = j;
				j3 = j+1;
				ksi1 = Globals::getY(j1);
				ksi2 = Globals::getY(j2);
				ksi3 = Globals::getY(j3);
			} 
			
			else {
				i1 = i2 = i3 = i-1;
				j1 = j;
				j2 = j+1;
				j3 = j+2;
				ksi1 = Globals::getY(j1);
				ksi2 = Globals::getY(j2);
				ksi3 = Globals::getY(j3);
			}
		}
	}
	
	// Allocate space for the upwind Intensity:
	rays[g][m].Iu = new (nothrow) double [Globals::getN()];
	
	// And now after all has been done, we interpolate:
	if (dointerpol == 1){
		rays[g][m].Su = interpolate_bezier (ksi1, ksi2, ksi3, Grid[i1][j1].S, Grid[i2][j2].S , Grid[i3][j3].S, ksi);
		rays[g][m].Chi = interpolate_bezier (ksi1, ksi2, ksi3, Grid[i1][j1].Chi, Grid[i2][j2].Chi , Grid[i3][j3].Chi, ksi);
		rays[g][m].Chi = 1.0;
		for (int n=0; n<Globals::getN(); n++)
			rays[g][m].Iu[n] = interpolate_bezier (ksi1, ksi2, ksi3, Grid[i1][j1].I[g][m][n], Grid[i2][j2].I[g][m][n], Grid[i3][j3].I[g][m][n], ksi);
	}
	
	else{

		if (periodic_boundary == 0){
			for (int n = 0; n<Globals::getN(); n++){
			
			// We now have four input directions:
			
				if (m<Globals::getNP_C(g)/4){
				
					rays[g][m].Iu[n] = I_inc_1;
				}
			
				else if (m<Globals::getNP_C(g)/2){
				
					rays[g][m].Iu[n] = I_inc_2;
				}
			
				else if (m<Globals::getNP_C(g) * 3/4){
				
					rays[g][m].Iu[n] = I_inc_3;
				}
			
				else if (m<Globals::getNP_C(g)){
				
					rays[g][m].Iu[n] = I_inc_4;
				}
			}
		}		
	}

	// Before we go out, we have to correct in the case the boundary conditions are periodic

	int n;
	int N = Globals::getN();
	//int NX = Globals::getNX();
	//int NY = Globals::getNY();
	double delta;
	double index_i;

	if (periodic_boundary == 1){

		if (j==0 && Globals::getphi_C(g,m) <= pi)
			for (n=0; n<N; n++)
				rays[g][m].Iu[n] = 0.0; // Incoming illumination equal to zero! 

		else if (j==Globals::getNY()-1 && Globals::getphi_C(g,m) >= pi)
			for (n=0; n<N; n++)
				rays[g][m].Iu[n] = B; // Incoming illumination equal to the local Planck function

		else if (i==0 && Globals::getcosphi_C(g,m) >= 0){ // Incoming illumation from the left is identical to the incoming illumination
				// from the left in the case when i == NX-1

			Grid[NX-1][j].allocate_rays();
			Grid[NX-1][j].compute_upwind_points(g, m, Grid);
			Grid[NX-1][j].interpolate_simple(g,m, Grid);

			for (n=0; n<N; n++)
				Grid[i][j].rays[g][m].Iu[n] = Grid[NX-1][j].rays[g][m].Iu[n];

			delete []Grid[NX-1][j].rays[g][m].Iu;

			Grid[NX-1][j].deallocate_rays();

		}

		else if (i==NX-1 && Globals::getcosphi_C(g,m) <= 0){ // Incoming illumation from the right is identical to the incoming illumination
				// from the right in the case when i == 0

			Grid[0][j].allocate_rays();
			Grid[0][j].compute_upwind_points(g, m, Grid);
			Grid[0][j].interpolate_simple(g,m, Grid);

			for (n=0; n<N; n++)
				Grid[i][j].rays[g][m].Iu[n] = Grid[0][j].rays[g][m].Iu[n];

			delete []Grid[0][j].rays[g][m].Iu;

			Grid[0][j].deallocate_rays();

		} 

	}

	if (periodic_boundary == 2){ // Simplest way to account for periodic boundaries. Very slow. Do not use.

		if (j==0 && Globals::getphi_C(g,m) <= pi)
			for (n=0; n<N; n++)
				rays[g][m].Iu[n] = 0.0; // Incoming illumination equal to zero! 

		else if (j==Globals::getNY()-1 && Globals::getphi_C(g,m) >= pi)
			for (n=0; n<N; n++)
				rays[g][m].Iu[n] = B; // Incoming illumination equal to the local Planck function 

		else if (i==0 && (Globals::getphi_C(g,m) <= pi/2 || Globals::getphi_C(g,m) >= 3*pi/2)){
			
			for (n=0; n<N; n++)
				rays[g][m].Iu[n] = Grid[Globals::getNX() -1][j].I[g][m][n]; // It is not terrible but is a bit slow.

		}

		else if (i==Globals::getNX()-1 && Globals::getphi_C(g,m) <= 3 * pi/2 && Globals::getphi_C(g,m) >= pi/2){
			for (n=0; n<N; n++)
				rays[g][m].Iu[n] = Grid[0][j].I[g][m][n]; // It is not terrible but is a bit slow.
		}

	}
	
	return 0;
}


// Computation of transport coefficients:

int Point::compute_pqr (int g, int m){
	
	// This is a function which computes p,q,r coefficients with the assumption that Spy and Spz are symmetrically computed 
	// according to the last direction of sweeping the grid
	
	double margin = 1e-4;
	
	int N = Globals::getN();
	rays[g][m].dt_u = new (nothrow) double [N];
	
	int n;
	
	// Now distance in 3D space:
	double distance = sqrt ((X - rays[g][m].Xu) * (X - rays[g][m].Xu) + (Y - rays[g][m].Yu) * (Y - rays[g][m].Yu)) / Globals::getcostheta(g);
	double dt;
	
	for (n=0; n<N; n++)
		rays[g][m].dt_u[n] = distance * Globals::getprofile(n);
		
	// Here we compute p,q,r with the assumption of parabolic shape of the source function:
		
	for (n=0; n<N; n++){
		
		dt = rays[g][m].dt_u[n];
		
		if (dt < e-20){
			p[g][m][n] = 0.0;
			q[g][m][n] = 0.0;
			r[g][m][n] = 0.0;
			exxp[g][m][n] = 1.0;
		}
		
		else if (dt < margin){
			
			p[g][m][n] = 2.0/3.0 * dt  - dt*dt/4.0 + dt*dt*dt/15.0;
			q[g][m][n] = dt / 3.0 - dt * dt / 4.0 + dt * dt * dt / 10.0;
			r[g][m][n] = -dt * dt / 6.0 + dt * dt * dt / 12.0 - dt * dt * dt * dt/40;
			exxp[g][m][n] = 1.0 - dt + dt * dt * 0.5 - dt * dt * dt / 6.0;
		}
		
		else {
			
			p[g][m][n] = 1.0 + 2.0/dt * exp(-dt) - 2.0/dt/dt * (1.0 - exp(-dt));
			q[g][m][n] = -exp(-dt) - 2.0/dt * exp(-dt) + 2.0/dt/dt * (1.0 - exp(-dt));
			r[g][m][n] = -1.0 - exp(-dt) + 2.0/dt * (1.0 - exp(-dt));
			exxp[g][m][n] = exp(-dt);
		}
	}
// ------------------------------------------------------------------------------------------------------------------------------

	// But now some modifications have to be made:
	
	// Then	 we modify for the contribution of the local source function to the x derivative:
	
	for (n=0; n<N; n++){
			
		if (i==0){
		
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / (Globals::getX(0) - Globals::getX(1));
			
		
		}
		 
		else if (i == Globals::getNX() -1){
	
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / (Globals::getX(i) - Globals::getX(i-1));
			
		}
			
		else{
			
			double x1 = Globals::getX(i-1);
			double x2 = Globals::getX(i);
			double x3 = Globals::getX(i+1);
		
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n)
				* 1.0 / (x3 - x1) * ((x3 - x2) / (x2 - x1) + (x1 - x2) / (x3 - x2));
			
		}
	}
	
		
// ---------------------------------------------------------------------------------------------------------------------------------
		
	// Then we modify for the contribution of the local source function to the y derivative:
	
	
	for (n=0; n<N; n++){
			
		if (j == 0){
		
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / ( Globals::getY(j) - Globals::getY(j+1) );
		
		}
			
		else if (j == Globals::getNY()-1){
		
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / ( Globals::getY(j) - Globals::getY(j-1) );
		
		}
			
		else {
			
			double y1 = Globals::getY(j-1);
			double y2 = Globals::getY(j);
			double y3 = Globals::getY(j+1);
			
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) *
				1.0 / (y3 - y1) * ((y3 - y2) / (y2 - y1) + (y1 - y2) / (y3 - y2));
				
		}	
	}
	
	// And compute rx, and ry:
		
	for (n=0; n<N; n++){
			
		rx[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n);
		ry[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n);
			
	}
	
	
	delete []rays[g][m].dt_u;
	
	return 0;
}

int Point::compute_pqr_4th_direction (int g, int m){
	
	// This is a function which computes p,q,r coefficients with the assumption that Spy and Spz are symmetrically computed 
	// according to the last direction of sweeping the grid
	
	double margin = 1e-4;
	
	int N = Globals::getN();
	rays[g][m].dt_u = new (nothrow) double [N];
	
	int n;
	
	// Now distance in 3D space:
	double distance = sqrt ((X - rays[g][m].Xu) * (X - rays[g][m].Xu) + (Y - rays[g][m].Yu) * (Y - rays[g][m].Yu)) / Globals::getcostheta(g);
	double dt;
	
	for (n=0; n<N; n++)
		rays[g][m].dt_u[n] = distance * Globals::getprofile(n);
		
	// Here we compute p,q,r with the assumption of parabolic shape of the source function:
		
	for (n=0; n<N; n++){
		
		dt = rays[g][m].dt_u[n];
		
		if (dt < e-20){
			p[g][m][n] = 0.0;
			q[g][m][n] = 0.0;
			r[g][m][n] = 0.0;
			exxp[g][m][n] = 1.0;
		}
		
		else if (dt < margin){
			
			p[g][m][n] = 2.0/3.0 * dt  - dt*dt/4.0 + dt*dt*dt/15.0;
			q[g][m][n] = dt / 3.0 - dt * dt / 4.0 + dt * dt * dt / 10.0;
			r[g][m][n] = -dt * dt / 6.0 + dt * dt * dt / 12.0 - dt * dt * dt * dt/40;
			exxp[g][m][n] = 1.0 - dt + dt * dt * 0.5 - dt * dt * dt / 6.0;
		}
		
		else {
			
			p[g][m][n] = 1.0 + 2.0/dt * exp(-dt) - 2.0/dt/dt * (1.0 - exp(-dt));
			q[g][m][n] = -exp(-dt) - 2.0/dt * exp(-dt) + 2.0/dt/dt * (1.0 - exp(-dt));
			r[g][m][n] = -1.0 - exp(-dt) + 2.0/dt * (1.0 - exp(-dt));
			exxp[g][m][n] = exp(-dt);
		}
	}
// ------------------------------------------------------------------------------------------------------------------------------

	// But now some modifications have to be made:
	
	// Then	 we modify for the contribution of the local source function to the x derivative:
	
	for (n=0; n<N; n++){
			
		if (i==0){
		
			p[g][m][n] += 0.0;
			
		}
		 
		else if (i == 1){
	
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / (Globals::getX(i) - Globals::getX(i-1));
			
		}
			
		else{
			
			double x1 = Globals::getX(i-2);
			double x2 = Globals::getX(i-1);
			double x3 = Globals::getX(i);
		
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n)
				* (1.0 / (x3 - x1) + 1.0 / (x3 - x2));
			
		}
	}
	
		
// ---------------------------------------------------------------------------------------------------------------------------------
		
	// Then we modify for the contribution of the local source function to the y derivative:
	
	
	for (n=0; n<N; n++){
			
		if (j == Globals::getNY()-1){
		
			p[g][m][n] += 0.0;
		
		}
			
		else if (j == Globals::getNY()-2){
		
			p[g][m][n] += - r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / ( Globals::getY(j+1) - Globals::getY(j) );
		
		}
			
		else {
			
			double y1 = Globals::getY(j);
			double y2 = Globals::getY(j+1);
			double y3 = Globals::getY(j+2);
			
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) *
				(1.0 / (y1 - y2) + 1.0 / (y1 - y3));
				
		}	
	}
	
	// And compute rx, and ry:
		
	for (n=0; n<N; n++){
			
		rx[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n);
		ry[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n);
			
	}
	
	
	delete []rays[g][m].dt_u;
	
	return 0;
}

int Point::compute_pqr_explicit_derivatives (int g, int m){
	
	// This is a function which computes p,q,r coefficients with the assumption that Spy and Spz are symmetrically computed 
	// according to the last direction of sweeping the grid
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	double margin = 1e-4;
	
	int N = Globals::getN();
	rays[g][m].dt_u = new (nothrow) double [N];
	
	int n;
	
	// Now distance in 3D space:
	double distance = sqrt ((X - rays[g][m].Xu) * (X - rays[g][m].Xu) + (Y - rays[g][m].Yu) * (Y - rays[g][m].Yu)) / Globals::getcostheta(g);
	double dt;
	
	for (n=0; n<N; n++)
		rays[g][m].dt_u[n] = distance * Globals::getprofile(n);
		
	// Here we compute p,q,r with the assumption of parabolic shape of the source function:
		
	for (n=0; n<N; n++){
		
		dt = rays[g][m].dt_u[n];
		
		if (dt < e-20){
			p[g][m][n] = 0.0;
			q[g][m][n] = 0.0;
			r[g][m][n] = 0.0;
			exxp[g][m][n] = 1.0;
		}
		
		else if (dt < margin){
			
			p[g][m][n] = 2.0/3.0 * dt  - dt*dt/4.0 + dt*dt*dt/15.0;
			q[g][m][n] = dt / 3.0 - dt * dt / 4.0 + dt * dt * dt / 10.0;
			r[g][m][n] = -dt * dt / 6.0 + dt * dt * dt / 12.0 - dt * dt * dt * dt/40;
			exxp[g][m][n] = 1.0 - dt + dt * dt * 0.5 - dt * dt * dt / 6.0;
		}
		
		else {
			
			p[g][m][n] = 1.0 + 2.0/dt * exp(-dt) - 2.0/dt/dt * (1.0 - exp(-dt));
			q[g][m][n] = -exp(-dt) - 2.0/dt * exp(-dt) + 2.0/dt/dt * (1.0 - exp(-dt));
			r[g][m][n] = -1.0 - exp(-dt) + 2.0/dt * (1.0 - exp(-dt));
			exxp[g][m][n] = exp(-dt);
		}
	}
// ------------------------------------------------------------------------------------------------------------------------------

	// But now some modifications have to be made:
	
	// Then	 we modify for the contribution of the local source function to the x derivative:
	
	for (n=0; n<N; n++){
			
		if (i==0){
		
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / (Globals::getX(0) - Globals::getX(1));
			
		
		}
		 
		else if (i == Globals::getNX() -1){
	
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / (Globals::getX(i) - Globals::getX(i-1));
			
		}
			
		else{
			
			double x1 = Globals::getX(i-1);
			double x2 = Globals::getX(i);
			double x3 = Globals::getX(i+1);
		
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n)
				* 1.0 / (x3 - x1) * ((x3 - x2) / (x2 - x1) + (x1 - x2) / (x3 - x2));
			
		}
	}
	
		
// ---------------------------------------------------------------------------------------------------------------------------------
		
	// Then we modify for the contribution of the local source function to the y derivative:
	
	
	for (n=0; n<N; n++){
			
		if (j == 0){
		
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / ( Globals::getY(j) - Globals::getY(j+1) );
		
		}
			
		else if (j == Globals::getNY()-1){
		
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / ( Globals::getY(j) - Globals::getY(j-1) );
		
		}
			
		else {
			
			double y1 = Globals::getY(j-1);
			double y2 = Globals::getY(j);
			double y3 = Globals::getY(j+1);
			
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) *
				1.0 / (y3 - y1) * ((y3 - y2) / (y2 - y1) + (y1 - y2) / (y3 - y2));
				
		}	
	}
	
	// And compute r1 .. r4 :
	
	double x1, x2, x3;
	double y1, y2, y3;
		
	for (n=0; n<N; n++){
		
		// First r1 and r2, which arise due to the derivative with respect to x:
		
		if (i == 0){
			x1 = Globals::getX(i);
			x2 = Globals::getX(i+1);
			r2[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) / (x2 - x1);
		}
			
		else if (i == NX-1){
			x1 = Globals::getX(i-1);
			x2 = Globals::getX(i);
			r1[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) / -(x2 - x1);
		}
			
		else {
			x1 = Globals::getX(i-1);
			x2 = Globals::getX(i);
			x3 = Globals::getX(i+1);
			r1[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) *
			(x2 - x3) / (x1 - x2) / (x1 - x3);
			r2[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) *
			(x2 - x1) / (x3 - x1) / (x3 - x2);
		} 
		
		// And then r3 and r4, which arise due to the derivative with respect to y:
		
		if (j == 0){
			y1 = Globals::getY(j);
			y2 = Globals::getY(j+1);
			r4[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) / (y2 - y1);
		}
			
		else if (j == NX-1){
			y1 = Globals::getY(j);
			y2 = Globals::getY(j+1);
			r3[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) / -(y2 - y1);
		}
		
		else {
			y1 = Globals::getY(j-1);
			y2 = Globals::getY(j);
			y3 = Globals::getY(j+1);
			r3[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) *
			(y2 - y3) / (y1 - y2) / (y1 - y3);
			r4[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) *
			(y2 - y1) / (y3 - y1) / (y3 - y2);
		}	
	}
	
	delete []rays[g][m].dt_u;
	
	return 0;
}

int Point::compute_pqr_explicit_full (int g, int m){
	
	// And ultimately this is the FULLY explicit version:
	
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	double margin = 0.01;
	
	int N = Globals::getN();
	rays[g][m].dt_u = new (nothrow) double [N];
	
	int n;
	
	// Now distance in 3D space:
	double distance = sqrt ((X - rays[g][m].Xu) * (X - rays[g][m].Xu) + (Y - rays[g][m].Yu) * (Y - rays[g][m].Yu)) / Globals::getcostheta(g);
	double dt=0.0;
	
	for (n=0; n<N; n++)
		rays[g][m].dt_u[n] = distance * (1.0 + Chi)* Globals::getprofile(n);
		
	// Here we compute p,q,r with the assumption of parabolic shape of the source function:
		
	for (n=0; n<N; n++){
		
		dt = rays[g][m].dt_u[n];
		
		if (dt < e-20){
			p[g][m][n] = 0.0;
			q[g][m][n] = 0.0;
			r[g][m][n] = 0.0;
			exxp[g][m][n] = 1.0;
		}
		
		else if (dt < margin){
			
			p[g][m][n] = 2.0/3.0 * dt  - dt*dt/4.0 + dt*dt*dt/15.0;
			q[g][m][n] = dt / 3.0 - dt * dt / 4.0 + dt * dt * dt / 10.0;
			r[g][m][n] = -dt * dt / 6.0 + dt * dt * dt / 12.0 - dt * dt * dt * dt/40;
			exxp[g][m][n] = 1.0 - dt + dt * dt * 0.5 - dt * dt * dt / 6.0;
		}
		
		else {
			
			p[g][m][n] = 1.0 + 2.0/dt * exp(-dt) - 2.0/dt/dt * (1.0 - exp(-dt));
			q[g][m][n] = -exp(-dt) - 2.0/dt * exp(-dt) + 2.0/dt/dt * (1.0 - exp(-dt));
			r[g][m][n] = -1.0 - exp(-dt) + 2.0/dt * (1.0 - exp(-dt));
			exxp[g][m][n] = exp(-dt);
		}
	}
// ------------------------------------------------------------------------------------------------------------------------------

	// But now some modifications have to be made:
	
	// Then	 we modify for the contribution of the local source function to the x derivative:

	double x1, x2, x3;
	
	for (n=0; n<N; n++){
			
		if (i==0){
		
			if (periodic_boundary == 1){

				x1 = Globals::getX(NX-2) - Globals::getX(NX-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(i+1);

				p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n)
					* 1.0 / (x3 - x1) * ((x3 - x2) / (x2 - x1) + (x1 - x2) / (x3 - x2)); 

			}
			
			else 
			
				p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / (Globals::getX(0) - Globals::getX(1));
			
		}
		 
		else if (i == Globals::getNX() -1){

			if (periodic_boundary == 1){

				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(NX-1) + Globals::getX(1);

				p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n)
					* 1.0 / (x3 - x1) * ((x3 - x2) / (x2 - x1) + (x1 - x2) / (x3 - x2));
				
			}

			else 
	
				p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / (Globals::getX(i) - Globals::getX(i-1));
			
		}
			
		else{
			
			x1 = Globals::getX(i-1);
			x2 = Globals::getX(i);
			x3 = Globals::getX(i+1);
		
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n)
				* 1.0 / (x3 - x1) * ((x3 - x2) / (x2 - x1) + (x1 - x2) / (x3 - x2));
			
		}
	}
	
		
// ---------------------------------------------------------------------------------------------------------------------------------
		
	// Then we modify for the contribution of the local source function to the y derivative:
	
	
	for (n=0; n<N; n++){
			
		if (j == 0){
		
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / ( Globals::getY(j) - Globals::getY(j+1) );
		
		}
			
		else if (j == Globals::getNY()-1){
		
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / ( Globals::getY(j) - Globals::getY(j-1) );
		
		}
			
		else {
			
			double y1 = Globals::getY(j-1);
			double y2 = Globals::getY(j);
			double y3 = Globals::getY(j+1);
			
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) *
				1.0 / (y3 - y1) * ((y3 - y2) / (y2 - y1) + (y1 - y2) / (y3 - y2));
				
		}	
	}
	
	// And compute r1 .. r4 :
	
	//double x1, x2, x3;
	double y1, y2, y3;
	double x, y;
		
	for (n=0; n<N; n++){
		
		// First r1 and r2, which arise due to the derivative with respect to x:
		
		if (i == 0){

			if (periodic_boundary == 1){

				x1 = Globals::getX(NX-2) - Globals::getX(NX-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(i+1);

				r1[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) *
					(x2 - x3) / (x1 - x2) / (x1 - x3);
				r2[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) *
					(x2 - x1) / (x3 - x1) / (x3 - x2);

			}

			else {
				x1 = Globals::getX(i);
				x2 = Globals::getX(i+1);
				r2[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) / (x2 - x1);
			}
		}
			
		else if (i == NX-1){

			if (periodic_boundary == 1){

				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(NX-1) + Globals::getX(1);

				r1[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) *
					(x2 - x3) / (x1 - x2) / (x1 - x3);
				r2[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) *
					(x2 - x1) / (x3 - x1) / (x3 - x2);

			}

			else {
				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				r1[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) / -(x2 - x1);
			}
		}
			
		else {
			x1 = Globals::getX(i-1);
			x2 = Globals::getX(i);
			x3 = Globals::getX(i+1);
			r1[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) *
			(x2 - x3) / (x1 - x2) / (x1 - x3);
			r2[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) *
			(x2 - x1) / (x3 - x1) / (x3 - x2);
		} 

		// And then r3 and r4, which arise due to the derivative with respect to y:
		
		if (j == 0){
			y1 = Globals::getY(j);
			y2 = Globals::getY(j+1);
			r4[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) / (y2 - y1);
		}
			
		else if (j == NY-1){
			y1 = Globals::getY(j);
			y2 = Globals::getY(j+1);
			r3[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) / -(y2 - y1);
		}
		
		else {
			y1 = Globals::getY(j-1);
			y2 = Globals::getY(j);
			y3 = Globals::getY(j+1);
			r3[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) *
			(y2 - y3) / (y1 - y2) / (y1 - y3);
			r4[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) *
			(y2 - y1) / (y3 - y1) / (y3 - y2);
		}	
	}
	
	// ----------------------------------------------------------------------------------------------------------------------------
	
	// And then modify them further due to the contribution of surrounding source functions to the upwind source function:
	
	// For turning on and off factors and contributions:
	
	double factor1, factor2, factor3, factor4;
	
	// 1: there is factor, 0 there is no factor.
	
	factor1 = 0.0;
	factor2 = 0.0;
	factor3 = 0.0;
	factor4 = 0.0;
	
	
	// QUADRANT 1:
	
	
	if (Globals::getphi_C(g, m) <= pi/2) {
		if (rays[g][m].type_u == 1){
			
			if (i == 0){

				if (periodic_boundary == 1){

					x1 = Globals::getX(NX-2) - Globals::getX(NX-1);
					x2 = Globals::getX(i);
					x3 = Globals::getX(i+1);
					x = rays[g][m].Xu;

					for (n=0; n<N; n++){
						r5[g][m][n] += (1.0 - factor1) * (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
						r3[g][m][n] += (1.0 - factor1) * (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
						r6[g][m][n] += (1.0 - factor1) * (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
					}

				}

				else {	
					for (n=0; n<N; n++){
						r5[g][m][n] += 0.0;
						r3[g][m][n] += 0.0;
						r6[g][m][n] += 0.0;
					}
				}
			}
			
			else if (i == NX-1){
				
				if (periodic_boundary == 1){

					x1 = Globals::getX(i-1);
					x2 = Globals::getX(i);
					x3 = Globals::getX(NX-1) + Globals::getX(1);	
					x = rays[g][m].Xu;

					for (n=0; n<N; n++){
						r5[g][m][n] += (1.0 - factor1) * (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
						r3[g][m][n] += (1.0 - factor1) * (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
						r6[g][m][n] += (1.0 - factor1) * (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
					}

				}

				else {
					x1 = Globals::getX(i-1);
					x2 = Globals::getX(i);
					x = rays[g][m].Xu;
				
					for (n=0; n<N; n++){
						r5[g][m][n] += (1.0 - factor1) * (x2 - x) / (x2 - x1) * q[g][m][n];
						r3[g][m][n] += (1.0 - factor1) * (x1 - x) / (x1 - x2) * q[g][m][n];
						r6[g][m][n] += 0.0;
					}
				}
			}
			
			else {
				
				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(i+1);
				x = rays[g][m].Xu;

				
				for (n=0; n<N; n++){
					r5[g][m][n] += (1.0 - factor1) * (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
					r3[g][m][n] += (1.0 - factor1) * (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
					r6[g][m][n] += (1.0 - factor1) * (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
				}
			}
		}
		
		else if (rays[g][m].type_u == 2){
			
			if (j == 0){
				for (n=0; n<N; n++){
					r5[g][m][n] += 0.0;
					r1[g][m][n] += 0.0;
					r7[g][m][n] += 0.0;
				}
			}
			
			else if (j == NY-1){
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += (1.0 - factor1) * (y2 - y) / (y2 - y1) * q[g][m][n];
					r1[g][m][n] += (1.0 - factor1) * (y1 - y) / (y1 - y2) * q[g][m][n];
					r7[g][m][n] += 0.0;
				}
			}
			
			else {
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y3 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += (1.0 - factor1) * (y2 - y) * (y3 - y) / (y2 - y1) / (y3 - y1) * q[g][m][n];
					r1[g][m][n] += (1.0 - factor1) * (y1 - y) * (y3 - y) / (y1 - y2) / (y3 - y2) * q[g][m][n];
					r7[g][m][n] += (1.0 - factor1) * (y1 - y) * (y2 - y) / (y1 - y3) / (y2 - y3) * q[g][m][n];
				}
			}
		}
	}
	
	// QUADRANT 2:
	
	else if (Globals::getphi_C(g, m) <= pi){
		if (rays[g][m].type_u == 1){
			
			if (i == NX - 1){

				if (periodic_boundary == 1){

					x1 = Globals::getX(i-1);
					x2 = Globals::getX(i);
					x3 = Globals::getX(NX-1) + Globals::getX(1);
					x = rays[g][m].Xu;
					
					for (n=0; n<N; n++){
						r5[g][m][n] += (1.0 - factor1) * (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
						r3[g][m][n] += (1.0 - factor1) * (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
						r6[g][m][n] += (1.0 - factor1) * (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
					}

				}

				else {
					for (n=0; n<N; n++){
						r5[g][m][n] += 0.0;
						r3[g][m][n] += 0.0;
						r6[g][m][n] += 0.0;
					}
				}
			}
			
			else if (i == 0){

				if (periodic_boundary == 1){

					x1 = Globals::getX(NX-2) - Globals::getX(NX-1);
					x2 = Globals::getX(i);
					x3 = Globals::getX(i+1);
					x = rays[g][m].Xu;

					for (n=0; n<N; n++){
						r5[g][m][n] += (1.0 - factor1) * (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
						r3[g][m][n] += (1.0 - factor1) * (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
						r6[g][m][n] += (1.0 - factor1) * (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
					}

				}

				else {
				
					x1 = Globals::getX(i);
					x2 = Globals::getX(i+1);
					x = rays[g][m].Xu;
				
					for (n=0; n<N; n++){
						r5[g][m][n] += 0.0;
						r3[g][m][n] += (1.0 - factor2) * (x2 - x) / (x2 - x1) * q[g][m][n];
						r6[g][m][n] += (1.0 - factor2) * (x1 - x) / (x1 - x2) * q[g][m][n];
					}
				}
			}
			
			else {
				
				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(i+1);
				x = rays[g][m].Xu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += (1.0 - factor2) * (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
					r3[g][m][n] += (1.0 - factor2) * (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
					r6[g][m][n] += (1.0 - factor2) * (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
				}
			}
		}
		
		else if (rays[g][m].type_u == 2){
			
			if (j == 0){
				for (n=0; n<N; n++){
					r6[g][m][n] += 0.0;
					r2[g][m][n] += 0.0;
					r8[g][m][n] += 0.0;
				}
			}
			
			else if (j == NY-1){
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r6[g][m][n] += (1.0 - factor2) * (y2 - y) / (y2 - y1) * q[g][m][n];
					r2[g][m][n] += (1.0 - factor2) * (y1 - y) / (y1 - y2) * q[g][m][n];
					r8[g][m][n] += 0.0;
				}
			}
			
			else {
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y3 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r6[g][m][n] += (1.0 - factor2) * (y2 - y) * (y3 - y) / (y2 - y1) / (y3 - y1) * q[g][m][n];
					r2[g][m][n] += (1.0 - factor2) * (y1 - y) * (y3 - y) / (y1 - y2) / (y3 - y2) * q[g][m][n];
					r8[g][m][n] += (1.0 - factor2) * (y1 - y) * (y2 - y) / (y1 - y3) / (y2 - y3) * q[g][m][n];
				}
			}
		}
	}
	
	// QUADRANT 3:
	
	else if (Globals::getphi_C(g, m) <= pi * 3 / 2){
		if (rays[g][m].type_u == 1){
			
			if (i == NX - 1){

				if (periodic_boundary == 1){

					x1 = Globals::getX(i-1);
					x2 = Globals::getX(i);
					x3 = Globals::getX(NX-1) + Globals::getX(1);
					x = rays[g][m].Xu;

					for (n=0; n<N; n++){
						r7[g][m][n] += (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
						r4[g][m][n] += (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
						r8[g][m][n] += (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
					}
				}

				else {
					for (n=0; n<N; n++){
						r7[g][m][n] += 0.0;
						r4[g][m][n] += 0.0;
						r8[g][m][n] += 0.0;
					}
				}
			}
			
			else if (i == 0 ){

				if (periodic_boundary == 1){

					x1 = Globals::getX(NX-2) - Globals::getX(NX-1);
					x2 = Globals::getX(i);
					x3 = Globals::getX(i+1);
					x = rays[g][m].Xu;

					for (n=0; n<N; n++){
						r7[g][m][n] += (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
						r4[g][m][n] += (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
						r8[g][m][n] += (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
					}
				}
				
				else {
					x1 = Globals::getX(i);
					x2 = Globals::getX(i+1);
					x = rays[g][m].Xu;
				
					for (n=0; n<N; n++){
						r7[g][m][n] += 0.0;
						r4[g][m][n] += (x2 - x) / (x2 - x1) * q[g][m][n];
						r8[g][m][n] += (x1 - x) / (x1 - x2) * q[g][m][n];
					}
				}
			}
			
			else {
				
				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(i+1);
				x = rays[g][m].Xu;
				
				for (n=0; n<N; n++){
					r7[g][m][n] += (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
					r4[g][m][n] += (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
					r8[g][m][n] += (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
				}
			}
		}
		
		else if (rays[g][m].type_u == 2){
			
			if (j == NY-1){
				for (n=0; n<N; n++){
					r6[g][m][n] += 0.0;
					r2[g][m][n] += 0.0;
					r8[g][m][n] += 0.0;
				}
			}
			
			else if (j == 0){
				
				y1 = Globals::getY(j);
				y2 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r6[g][m][n] += 0.0;
					r2[g][m][n] += (y2 - y) / (y2 - y1) * q[g][m][n];
					r8[g][m][n] += (y1 - y) / (y1 - y2) * q[g][m][n];
				}
				
			}
			
			else {
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y3 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r6[g][m][n] += (y2 - y) * (y3 - y) / (y2 - y1) / (y3 - y1) * q[g][m][n];
					r2[g][m][n] += (y1 - y) * (y3 - y) / (y1 - y2) / (y3 - y2) * q[g][m][n];
					r8[g][m][n] += (y1 - y) * (y2 - y) / (y1 - y3) / (y2 - y3) * q[g][m][n];
				}
			}
		}
	}
	
	// QUADRANT 4:
	
	else if (Globals::getphi_C(g, m) <= pi * 4 / 2){
		if (rays[g][m].type_u == 1){
			
			if (i == 0){

				if (periodic_boundary == 1){

					x1 = Globals::getX(NX-2) - Globals::getX(NX-1);
					x2 = Globals::getX(i);
					x3 = Globals::getX(i+1);
					x = rays[g][m].Xu;

					for (n=0; n<N; n++){
						r7[g][m][n] += (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
						r4[g][m][n] += (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
						r8[g][m][n] += (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];

					}
				}
				

				else {
					for (n=0; n<N; n++){
						r7[g][m][n] += 0.0;
						r4[g][m][n] += 0.0;
						r8[g][m][n] += 0.0;
					}
				}
			}
			
			else if (i == NX-1){

				if (periodic_boundary == 1){

					x1 = Globals::getX(i-1);
					x2 = Globals::getX(i);
					x3 = Globals::getX(NX-1) + Globals::getX(1);
					x = rays[g][m].Xu;

					for (n=0; n<N; n++){
						r7[g][m][n] += (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
						r4[g][m][n] += (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
						r8[g][m][n] += (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];

					}

				}
				
				else {

					x1 = Globals::getX(i-1);
					x2 = Globals::getX(i);
					x = rays[g][m].Xu;
				
					for (n=0; n<N; n++){
						r7[g][m][n] += (x2 - x) / (x2 - x1) * q[g][m][n];
						r4[g][m][n] += (x1 - x) / (x1 - x2) * q[g][m][n];
						r8[g][m][n] += 0.0;
					}
				}
			}
			
			else {
				
				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(i+1);
				x = rays[g][m].Xu;
				
				for (n=0; n<N; n++){
					r7[g][m][n] += (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
					r4[g][m][n] += (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
					r8[g][m][n] += (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
				}
			}
		}
		
		else if (rays[g][m].type_u == 2){
			
			if (j == NY-1){
				for (n=0; n<N; n++){
					r5[g][m][n] += 0.0;
					r1[g][m][n] += 0.0;
					r7[g][m][n] += 0.0;
				}
			}
			
			else if (j == 0){
				
				y1 = Globals::getY(j);
				y2 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += 0.0;
					r1[g][m][n] += (y2 - y) / (y2 - y1) * q[g][m][n];
					r7[g][m][n] += (y1 - y) / (y1 - y2) * q[g][m][n];
				}
			}
			
			else {
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y3 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += (y2 - y) * (y3 - y) / (y2 - y1) / (y3 - y1) * q[g][m][n];
					r1[g][m][n] += (y1 - y) * (y3 - y) / (y1 - y2) / (y3 - y2) * q[g][m][n];
					r7[g][m][n] += (y1 - y) * (y2 - y) / (y1 - y3) / (y2 - y3) * q[g][m][n];
				}
			}
		}
	}
	
	delete []rays[g][m].dt_u;
	
	return 0;
}

// Formal solution related stuff:

int Point::formal_solution (int g, int m){
	
	
	// Here we compute the intensity:
	
	for (int n = 0; n<Globals::getN(); n++)
		I[g][m][n] = rays[g][m].Iu[n] * exxp[g][m][n] + p[g][m][n] * S + q[g][m][n] * rays[g][m].Su + rx[g][m][n] * Spx_mod + ry[g][m][n] * Spy_mod;
			
	// But in the heart of FBILI method is the computation of the AA and BB coefficients:
	
	// For the moment let's compute them Jacobi style:
	
	for (int n = 0; n<Globals::getN(); n++){
		
		AA += (rays[g][m].Iu[n] * exxp[g][m][n] + q[g][m][n] * rays[g][m].Su) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		BB += p[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		CX += rx[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		CY += ry[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
	}
	
	// Delete the upwind intensity as it is no longer needed:
	
	delete []rays[g][m].Iu;
	
	return 0;
}

int Point::formal_solution_explicit_derivatives (int g, int m, Point ** Grid, int iter){
	
	// First establish four neighboring source functions:
	
	double S1, S2, S3, S4;
	
	// First "x" ones:
	
	if (i==0){
		S1 = 0;
		S2 = Grid[i+1][j].S;
	}
	
	else if (i== Globals::getNX() - 1){
		S1 = Grid[i-1][j].S;
		S2 = 0;
	}
	
	else {
		S1 = Grid[i-1][j].S;
		S2 = Grid[i+1][j].S;
	}
	
	
	// Now "y" ones:
	
	if (j==0){
		S3 = 0;
		S4 = Grid[i][j+1].S;
	}
	
	else if (j == Globals::getNY() - 1){
		S3 = Grid[i][j-1].S;
		S4 = 0;
	}
	
	else {
		S3 = Grid[i][j-1].S;
		S4 = Grid[i][j+1].S;
	}
	
	// Compute the intensity and add appropriate coefficients:
	
	for (int n = 0; n<Globals::getN(); n++){
		
		I[g][m][n] = rays[g][m].Iu[n] * exxp[g][m][n] + p[g][m][n] * S + q[g][m][n] * rays[g][m].Su + 
			r1[g][m][n] * S1 + r2[g][m][n] * S2 + r3[g][m][n] * S3 + r4[g][m][n] * S4;
			
		//AA += (rays[g][m].Iu[n] * exxp[g][m][n] + q[g][m][n] * rays[g][m].Su) * 
			//0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			
		if (m < Globals::getNP_C(g) * 1 / 4)
			AA1 += (rays[g][m].Iu[n] * exxp[g][m][n] + q[g][m][n] * rays[g][m].Su) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		else if (m < Globals::getNP_C(g) * 2 / 4)
			AA2 += (rays[g][m].Iu[n] * exxp[g][m][n] + q[g][m][n] * rays[g][m].Su) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		else if (m < Globals::getNP_C(g) * 3 / 4)
			AA3 += (rays[g][m].Iu[n] * exxp[g][m][n] + q[g][m][n] * rays[g][m].Su) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		else if (m < Globals::getNP_C(g) * 4 / 4)
			AA4 += (rays[g][m].Iu[n] * exxp[g][m][n] + q[g][m][n] * rays[g][m].Su) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			 
		if (iter == 0){
		BB += p[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C1 += r1[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C2 += r2[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C3 += r3[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C4 += r4[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
	}
	
	delete []rays[g][m].Iu;
	
	return 0;
}

int Point::formal_solution_explicit_derivatives_factors (int g, int m, Point ** Grid, int iter){
	
	// First establish four neighboring source functions:
	
	double S1, S2, S3, S4;
	
	// First "x" ones:
	
	if (i==0){
		S1 = 0;
		S2 = Grid[i+1][j].S;
	}
	
	else if (i== Globals::getNX() - 1){
		S1 = Grid[i-1][j].S;
		S2 = 0;
	}
	
	else {
		S1 = Grid[i-1][j].S;
		S2 = Grid[i+1][j].S;
	}
	
	
	// Now "y" ones:
	
	if (j==0){
		S3 = 0;
		S4 = Grid[i][j+1].S;
	}
	
	else if (j == Globals::getNY() - 1){
		S3 = Grid[i][j-1].S;
		S4 = 0;
	}
	
	else {
		S3 = Grid[i][j-1].S;
		S4 = Grid[i][j+1].S;
	}
	
	// Compute intensity and add appropriate coefficients:
	
	for (int n = 0; n<Globals::getN(); n++){
		
		I[g][m][n] = rays[g][m].Iu[n] * exxp[g][m][n] + p[g][m][n] * S + q[g][m][n] * rays[g][m].Su + 
			r1[g][m][n] * S1 + r2[g][m][n] * S2 + r3[g][m][n] * S3 + r4[g][m][n] * S4;
				
		if (m < Globals::getNP_C(g) * 1 / 4){
			
			AA1 += (q[g][m][n] * rays[g][m].Su) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			
			BB1 += (rays[g][m].Iu[n] * exxp[g][m][n] / S + p[g][m][n]) *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 2 / 4){
			
			AA2 += (q[g][m][n] * rays[g][m].Su) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			
			BB2 += (rays[g][m].Iu[n] * exxp[g][m][n] / S + p[g][m][n]) *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 3 / 4){
			
			AA3 += (rays[g][m].Iu[n] * exxp[g][m][n] + q[g][m][n] * rays[g][m].Su) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			
			BB3 += p[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 4 / 4){
			
			AA4 += (rays[g][m].Iu[n] * exxp[g][m][n] + q[g][m][n] * rays[g][m].Su) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			
			BB4 += p[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
			 
		if (iter == 0){
		
		C1 += r1[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C2 += r2[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C3 += r3[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C4 += r4[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
	}
	
	delete []rays[g][m].Iu;
	
	return 0;
}

int Point::formal_solution_explicit_full (int g, int m, Point ** Grid, int iter){
	
		// First establish four neighboring source functions:
	
	double S1, S2, S3, S4, S5, S6, S7, S8;
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	// First "x" ones:
	
	if (i==0){
		S1 = 0;
		S2 = Grid[i+1][j].S;
	}
	
	else if (i== Globals::getNX() - 1){
		S1 = Grid[i-1][j].S;
		S2 = 0;
	}
	
	else {
		S1 = Grid[i-1][j].S;
		S2 = Grid[i+1][j].S;
	}
	
	
	// Now "y" ones:
	
	if (j==0){
		S3 = 0;
		S4 = Grid[i][j+1].S;
	}
	
	else if (j == Globals::getNY() - 1){
		S3 = Grid[i][j-1].S;
		S4 = 0;
	}
	
	else {
		S3 = Grid[i][j-1].S;
		S4 = Grid[i][j+1].S;
	}
	
	// And corner ones:
	
	if (i == 0 || j == 0)
		S5 = 0.0;
	
	else 
		S5 = Grid[i-1][j-1].S;

	if (i== 0 || j == NY-1)
		S7 = 0.0;
	else 
		S7 = Grid[i-1][j+1].S;
		
	if (i == NX - 1 || j == 0)
		S6 = 0.0;
	else 
		S6 = Grid[i+1][j-1].S;
		
	if (i == NX - 1 || j == NY-1)
		S8 = 0.0;
	else 
		S8 = Grid[i+1][j+1].S;

	////////////////////////////////////////////////////////////////////////////////////
	// A significant modifications in the case of periodic boundary conditions!!!!!!!!!!
	////////////////////////////////////////////////////////////////////////////////////

	if (periodic_boundary == 1){

		if (i == 0){ // Left boundary

			S1 = Grid[NX-2][j].S;
			
			if (j>0)
				S5 = Grid[NX-2][j-1].S;
			
			if (j<NY-1)
				S7 = Grid[NX-2][j+1].S;

		}

		if (i == Globals::getNX()-1){ // Right boundary

			S2 = Grid[1][j].S;

			if (j>0)
				S6 = Grid[1][j-1].S;

			if (j<NY-1)
				S8 = Grid[1][j+1].S;

		}
				
	}
	
	// Compute intensity and add appropriate coefficients:
	
	for (int n = 0; n<Globals::getN(); n++){
		
		I[g][m][n] = rays[g][m].Iu[n] * exxp[g][m][n] + p[g][m][n] * S  + 
			r1[g][m][n] * S1 + r2[g][m][n] * S2 + r3[g][m][n] * S3 + r4[g][m][n] * S4 + 
			r5[g][m][n] * S5 + r6[g][m][n] * S6 + r7[g][m][n] * S7 + r8[g][m][n] * S8;
			
		if (m < Globals::getNP_C(g) * 1 / 4){
			AA1 += (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB1 += p[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 2 / 4){
			AA2 += (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB2 += p[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 3 / 4){
			AA3 += (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB3 += p[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 4 / 4){
			AA4 += (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB4 += p[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
			 
		if (iter == 0){
		
		C1 += r1[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C2 += r2[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C3 += r3[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C4 += r4[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C5 += r5[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C6 += r6[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C7 += r7[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C8 += r8[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
	}
	
	delete []rays[g][m].Iu;
	
	return 0;
}

int Point::formal_solution_explicit_full_factors (int g, int m, Point ** Grid, int iter){
	
	// First establish four neighboring source functions:
	
	double S1, S2, S3, S4, S5, S6, S7, S8;
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	// First "x" ones:
	
	if (i==0){
		S1 = 0;
		S2 = Grid[i+1][j].S;
	}
	
	else if (i== Globals::getNX() - 1){
		S1 = Grid[i-1][j].S;
		S2 = 0;
	}
	
	else {
		S1 = Grid[i-1][j].S;
		S2 = Grid[i+1][j].S;
	}
	
	
	// Now "y" ones:
	
	if (j==0){
		S3 = 0;
		S4 = Grid[i][j+1].S;
	}
	
	else if (j == Globals::getNY() - 1){
		S3 = Grid[i][j-1].S;
		S4 = 0;
	}
	
	else {
		S3 = Grid[i][j-1].S;
		S4 = Grid[i][j+1].S;
	}
	
	// And corner ones:
	
	if (i == 0 || j == 0)
		S5 = 0.0;
	
	else 
		S5 = Grid[i-1][j-1].S;

	if (i== 0 || j == NY-1)
		S7 = 0.0;
	else 
		S7 = Grid[i-1][j+1].S;
		
	if (i == NX - 1 || j == 0)
		S6 = 0.0;
	else 
		S6 = Grid[i+1][j-1].S;
		
	if (i == NX - 1 || j == NY-1)
		S8 = 0.0;
	else 
		S8 = Grid[i+1][j+1].S;

	if (periodic_boundary == 1){

		if (i == 0){ // Left boundary

			S1 = Grid[NX-2][j].S;
			
			if (j>0)
				S5 = Grid[NX-2][j-1].S;
			
			if (j<NY-1)
				S7 = Grid[NX-2][j+1].S;

		}

		if (i == Globals::getNX()-1){ // Right boundary

			S2 = Grid[1][j].S;

			if (j>0)
				S6 = Grid[1][j-1].S;

			if (j<NY-1)
				S8 = Grid[1][j+1].S;

		}
			
	}
		
	// Compute intensity and add appropriate coefficients:
	
	for (int n = 0; n<Globals::getN(); n++){
		
		I[g][m][n] = rays[g][m].Iu[n] * exxp[g][m][n] + p[g][m][n] * S  + 
			r1[g][m][n] * S1 + r2[g][m][n] * S2 + r3[g][m][n] * S3 + r4[g][m][n] * S4 + 
			r5[g][m][n] * S5 + r6[g][m][n] * S6 + r7[g][m][n] * S7 + r8[g][m][n] * S8;
			
		if (m < Globals::getNP_C(g) * 1 / 4){
			AA1 += 0.0 * (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB1 += (p[g][m][n] + 1.0 * rays[g][m].Iu[n] * exxp[g][m][n] / S + 0.0 * q[g][m][n] * rays[g][m].Su / S) *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 2 / 4){
			AA2 += 1.0 * (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB2 += (p[g][m][n] + 0.0 * rays[g][m].Iu[n] * exxp[g][m][n] / S + 0.0 * q[g][m][n] * rays[g][m].Su / S) *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 3 / 4){
			AA3 += 1.0 * (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB3 += (p[g][m][n] + 0.0 * rays[g][m].Iu[n] * exxp[g][m][n] / S) *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 4 / 4){
			AA4 += (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB4 += p[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
			 
		if (iter == 0){
		
		C1 += r1[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C2 += r2[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C3 += r3[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C4 += r4[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C5 += r5[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C6 += r6[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C7 += r7[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C8 += r8[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
	}
	
	delete []rays[g][m].Iu;
	
	return 0;
}



// Intensity allocation/deallocation, and computation

int Point::allocate_I (){
	
	I = new (nothrow) double ** [Globals::getNT()];
	
	for (int g = 0; g<Globals::getNT(); g++){
		
		I[g] = new (nothrow) double * [Globals::getNP_C(g)];
		
		for (int m =0; m<Globals::getNP_C(g); m++){
			I[g][m] = new (nothrow) double [Globals::getN()];
			
			for (int n = 0; n<Globals::getN(); n++)
				I[g][m][n] = B;
			}
	}
	
	return 0;
}

int Point::deallocate_I (){
	
	for (int g = 0; g<Globals::getNT(); g++){
		
		for (int m = 0; m<Globals::getNP_C(g); m++)
			
			delete []I[g][m];
			
		delete []I[g];
	}
	
	delete []I;
	
	return 0;
}

// Derivatives:

int Point::compute_source_derivatives(Point ** Grid) {
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	double d1, d2, h1 ,h2;
	
	// First x derivative:
	
	if (i == 0)
		Spx = (Grid[1][j].S - S) / (Grid[1][j].X - X);
		
	else if (i == NX-1)
		Spx = (S - Grid[i-1][j].S) / (X - Grid[i-1][j].X);
		
	else {
		h1 = X - Grid[i-1][j].X;
		h2 = Grid[i+1][j].X - X;
		
		d1 = (S - Grid[i-1][j].S) / h1;
		d2 = (Grid[i+1][j].S - S) / h2;
		
		Spx = (h2 * d1 + h1 * d2) / (h1 + h2);
	} 
	
	// And then the y derivative:
	
	if (j == 0)
		Spy = (Grid[i][j+1].S - S) / (Grid[i][j+1].Y - Y);
		
	else if (j == NY-1)
		Spy = (S - Grid[i][j-1].S) / (Y - Grid[i][j-1].Y);
		
	else {
		h1 = Y - Grid[i][j-1].Y;
		h2 = Grid[i][j+1].Y - Y;
		
		d1 = (S - Grid[i][j-1].S) / h1;
		d2 = (Grid[i][j+1].S - S) / h2;
		
		Spy = (h2 * d1 + h1 * d2) / (h1 + h2);
	}

	return 0;
}

int Point::compute_source_derivatives_mod (Point ** Grid){
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	double d1, d2, h1 ,h2;
	
	// First x derivative:
	
	if (i == 0){
		Spx = (Grid[1][j].S - S) / (Grid[1][j].X - X);
		Spx_mod = (Grid[1][j].S) / (Grid[1][j].X - X);
	}
		
	else if (i == NX-1){
		Spx = (S - Grid[i-1][j].S) / (X - Grid[i-1][j].X);
		Spx = (- Grid[i-1][j].S) / (X - Grid[i-1][j].X);
	}
		
	else {
		h1 = X - Grid[i-1][j].X;
		h2 = Grid[i+1][j].X - X;
		
		d1 = (S - Grid[i-1][j].S) / h1;
		d2 = (Grid[i+1][j].S - S) / h2;
		
		Spx = (h2 * d1 + h1 * d2) / (h1 + h2);
		
		d1 = (- Grid[i-1][j].S) / h1;
		d2 = (Grid[i+1][j].S) / h2;
		Spx_mod = (h2 * d1 + h1 * d2) / (h1 + h2);
		
	} 
	
	// And then the y derivative:
	
	if (j == 0){
		Spy = (Grid[i][j+1].S - S) / (Grid[i][j+1].Y - Y);
		Spy_mod = (Grid[i][j+1].S) / (Grid[i][j+1].Y - Y);
	}
		
	else if (j == NY-1){
		Spy = (S - Grid[i][j-1].S) / (Y - Grid[i][j-1].Y);
		Spy_mod = (- Grid[i][j-1].S) / (Y - Grid[i][j-1].Y);
	}
		
	else {
		h1 = Y - Grid[i][j-1].Y;
		h2 = Grid[i][j+1].Y - Y;
		
		d1 = (S - Grid[i][j-1].S) / h1;
		d2 = (Grid[i][j+1].S - S) / h2;
		
		Spy = (h2 * d1 + h1 * d2) / (h1 + h2);
		
		d1 = (- Grid[i][j-1].S) / h1;
		d2 = (Grid[i][j+1].S) / h2;
		Spy_mod = (h2 * d1 + h1 * d2) / (h1 + h2);
	}
	
	return 0;
}

int Point::modify_derivatives_asymmetric (Point ** Grid){
	
	// This is the function which computes modified derivatives with asymmetric approximation:
	
	int NY = Globals::getNY();
	
	// First x:
	
	if (i == 0)
		Spx_mod = 0.0;
		
	else if (i == 1)
		Spx_mod = - Grid[i-1][j].S / (Globals::getX(i) - Globals::getX(i-1));
		
	else {
		
		double x1 = Globals::getX(i-2);
		double x2 = Globals::getX(i-1);
		double x3 = Globals::getX(i);
		
		Spx_mod = (x3 - x2) / (x1 - x2) / (x1 - x3) * Grid[i-2][j].S + (x3 - x1) / (x2 - x1) / (x2 - x3) * Grid[i-1][j].S;
	}
	
	// Then y:
	
	if (j == NY-1)
		Spy_mod = 0.0;
		
	else if (j == NY-2)
		Spy_mod = Grid[i][j+1].S / (Globals::getY(j+1) - Globals::getY(j));
		
	else {
		
		double y1 = Globals::getY(j);
		double y2 = Globals::getY(j+1);
		double y3 = Globals::getY(j+2);
		
		Spy_mod = (y1 - y3) / (y2 - y1) / (y2 - y3) * Grid[i][j+1].S + (y1 - y2) / (y3 - y1) / (y3 - y2) * Grid[i][j+2].S;
	}
	
	return 0;
}

// Source function computations:

int Point::compute_S_simple (){
	
	// Assuming that AA and BB are already computed we just compute dS and S and then reset AA and BB


	//if (i==0 && j==0){
	//	for (int g=0; g<Globals::getNT(); g++)
	//		for (int m=0; m<Globals::getNP_C(g); m++)
	//			for (int n=0; n<1; n++)
	//				cout << g << " " << m << " " << n << " " << I[g][m][n] << endl;
	//	}
	
	dS = (((1-ep) * AA + ep * B)) / (1 - (1-ep) * BB) - S;
	 
	S += dS;
	
	return 0;
}

int Point::compute_S_simple_explicit_derivatives (Point ** Grid){
	
	// First find neighboring source functions:
	
	double S1, S2, S3, S4;
	
	// First "x" ones:
	
	if (i==0){
		S1 = 0;
		S2 = Grid[i+1][j].S;
	}
	
	else if (i== Globals::getNX() - 1){
		S1 = Grid[i-1][j].S;
		S2 = 0;
	}
	
	else {
		S1 = Grid[i-1][j].S;
		S2 = Grid[i+1][j].S;
	}
	
	
	// Now "y" ones:
	
	if (j==0){
		S3 = 0;
		S4 = Grid[i][j+1].S;
	}
	
	else if (j == Globals::getNY() - 1){
		S3 = Grid[i][j-1].S;
		S4 = 0;
	}
	
	else {
		S3 = Grid[i][j-1].S;
		S4 = Grid[i][j+1].S;
	}
	
	AA = AA1 + AA2 + AA3 + AA4;
	
	AA += C1 * S1 + C2 * S2 + C3 * S3 + C4 * S4;
	
	AA *= 1.0;
	BB *= 1.0;
	
	dS = 1.0 * ((((1-ep) * AA + ep * B)) / (1 - (1-ep) * BB) - S);
	
	//if (i==32 && j == 0)
		//cout << "!!!!" << AA << " " << BB << " " << S << endl;
	 
	S += dS;
	
	return 0;
}

int Point::compute_S_simple_explicit_derivatives_factors (Point ** Grid){
		
	// First find neighboring source functions:
	
	double S1, S2, S3, S4;
	
	// First "x" ones:
	
	if (i==0){
		S1 = 0;
		S2 = Grid[i+1][j].S;
	}
	
	else if (i== Globals::getNX() - 1){
		S1 = Grid[i-1][j].S;
		S2 = 0;
	}
	
	else {
		S1 = Grid[i-1][j].S;
		S2 = Grid[i+1][j].S;
	}
	
	
	// Now "y" ones:
	
	if (j==0){
		S3 = 0;
		S4 = Grid[i][j+1].S;
	}
	
	else if (j == Globals::getNY() - 1){
		S3 = Grid[i][j-1].S;
		S4 = 0;
	}
	
	else {
		S3 = Grid[i][j-1].S;
		S4 = Grid[i][j+1].S;
	}
	
	AA = AA1 + AA2 + AA3 + AA4;
	BB = BB1 + BB2 + BB3 + BB4;
	
	AA += C1 * S1 + C2 * S2 + C3 * S3 + C4 * S4;
	
	AA *= 1.0;
	BB *= 1.0;
	
	dS = 1.0 * ((((1-ep) * AA + ep * B)) / (1 - (1-ep) * BB) - S);
	 
	S += dS;
	
	return 0;
}

int Point::compute_S_simple_explicit_full (Point ** Grid){
	
	// First find neighboring source functions:
	
	double S1, S2, S3, S4;
	double S5, S6, S7, S8;
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	// First "x" ones:
	
	if (i==0){
		S1 = 0;
		S2 = Grid[i+1][j].S;
	}
	
	else if (i== Globals::getNX() - 1){
		S1 = Grid[i-1][j].S;
		S2 = 0;
	}
	
	else {
		S1 = Grid[i-1][j].S;
		S2 = Grid[i+1][j].S;
	}
	
	
	// Now "y" ones:
	
	if (j==0){
		S3 = 0;
		S4 = Grid[i][j+1].S;
	}
	
	else if (j == Globals::getNY() - 1){
		S3 = Grid[i][j-1].S;
		S4 = 0;
	}
	
	else {
		S3 = Grid[i][j-1].S;
		S4 = Grid[i][j+1].S;
	}
	
	// And corner ones:
	
	if (i == 0 || j == 0)
		S5 = 0.0;
	
	else 
		S5 = Grid[i-1][j-1].S;

	if (i== 0 || j == NY-1)
		S7 = 0.0;
	else 
		S7 = Grid[i-1][j+1].S;
		
	if (i == NX - 1 || j == 0)
		S6 = 0.0;
	else 
		S6 = Grid[i+1][j-1].S;
		
	if (i == NX - 1 || j == NY-1)
		S8 = 0.0;
	else 
		S8 = Grid[i+1][j+1].S;


	if (periodic_boundary == 1){

		if (i == 0){ // Left boundary

			S1 = Grid[NX-2][j].S;
			
			if (j>0)
				S5 = Grid[NX-2][j-1].S;
			
			if (j<NY-1)
				S7 = Grid[NX-2][j+1].S;

		}

		if (i == Globals::getNX()-1){ // Right boundary

			S2 = Grid[1][j].S;

			if (j>0)
				S6 = Grid[1][j-1].S;

			if (j<NY-1)
				S8 = Grid[1][j+1].S;

		}
			
	}

	int factors = Globals::get_factor();
	double w = 0.8;

	if (factors == 1){
		BB1 += w * AA1/S;
		AA1 -= w * AA1;

		BB2 += w * AA2/S;
		AA2 -= w * AA2;

		//BB3 += w * AA3/S;
		//AA3 -= w * AA3;
	}
		
	
	AA = AA1 + AA2 + AA3 + AA4;
	
	AA += C1 * S1 + C2 * S2 + C3 * S3 + C4 * S4 + C5 * S5 + C6 * S6 + C7 * S7 + C8 * S8;
	BB = BB1 + BB2 + BB3 + BB4;
	
	dS = Globals::get_w() * ((((1-ep) * AA + ep * B)) / (1 - (1-ep) * BB) - S);
	
	if (S+dS < ep * B){
		dS = ep * B - S;
	} 
	 
	S += dS;
	
	return 0;
}

int Point::add_a_and_c(Point ** Grid){
	
	// First find neighboring source functions:
	
	double S1, S2, S3, S4;
	
	// First "x" ones:
	
	if (i==0){
		S1 = 0;
		S2 = Grid[i+1][j].S;
	}
	
	else if (i== Globals::getNX() - 1){
		S1 = Grid[i-1][j].S;
		S2 = 0;
	}
	
	else {
		S1 = Grid[i-1][j].S;
		S2 = Grid[i+1][j].S;
	}
	
	
	// Now "y" ones:
	
	if (j==0){
		S3 = 0;
		S4 = Grid[i][j+1].S;
	}
	
	else if (j == Globals::getNY() - 1){
		S3 = Grid[i][j-1].S;
		S4 = 0;
	}
	
	else {
		S3 = Grid[i][j-1].S;
		S4 = Grid[i][j+1].S;
	}

	
	AA = AA1 + AA2 + AA3 + AA4;
	
	AA += C1 * S1 + C2 * S2 + C3 * S3 + C4 * S4;

	BB = BB1 + BB2 + BB3 + BB4;
	
	return 0;
}

int Point::modify_a(Point ** Grid){

	double S1, S2, S3, S4;
	double S5, S6, S7, S8;
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	// First "x" ones:
	
	if (i==0){
		S1 = 0;
		S2 = Grid[i+1][j].S;
	}
	
	else if (i== Globals::getNX() - 1){
		S1 = Grid[i-1][j].S;
		S2 = 0;
	}
	
	else {
		S1 = Grid[i-1][j].S;
		S2 = Grid[i+1][j].S;
	}
	
	
	// Now "y" ones:
	
	if (j==0){
		S3 = 0;
		S4 = Grid[i][j+1].S;
	}
	
	else if (j == Globals::getNY() - 1){
		S3 = Grid[i][j-1].S;
		S4 = 0;
	}
	
	else {
		S3 = Grid[i][j-1].S;
		S4 = Grid[i][j+1].S;
	}
	
	// And corner ones:
	
	if (i == 0 || j == 0)
		S5 = 0.0;
	
	else 
		S5 = Grid[i-1][j-1].S;

	if (i== 0 || j == NY-1)
		S7 = 0.0;
	else 
		S7 = Grid[i-1][j+1].S;
		
	if (i == NX - 1 || j == 0)
		S6 = 0.0;
	else 
		S6 = Grid[i+1][j-1].S;
		
	if (i == NX - 1 || j == NY-1)
		S8 = 0.0;
	else 
		S8 = Grid[i+1][j+1].S;

	
	// And the modification for the case of the periodic boundary:

	if (periodic_boundary == 1){

		if (i == 0){ // Left boundary

			S1 = Grid[NX-2][j].S;
			
			if (j>0)
				S5 = Grid[NX-2][j-1].S;
			
			if (j<NY-1)
				S7 = Grid[NX-2][j+1].S;

		}

		if (i == Globals::getNX()-1){ // Right boundary

			S2 = Grid[1][j].S;

			if (j>0)
				S6 = Grid[1][j-1].S;

			if (j<NY-1)
				S8 = Grid[1][j+1].S;

		}
			
	}

	int factors = Globals::get_factor();
	double w = 1.0;

	if (factors == 1){
		BB1 += w * AA1/S;
		AA1 -= w * AA1;

		BB2 += w * AA2/S;
		AA2 -= w * AA2;
	}
		
	
	AA = AA1 + AA2 + AA3 + AA4;
	
	AA += C1 * S1 + C2 * S2 + C3 * S3 + C4 * S4 + C5 * S5 + C6 * S6 + C7 * S7 + C8 * S8;
	BB = BB1 + BB2 + BB3 + BB4;
	
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// EXPERIMENTAL FUNCTIONS ADDED AFTER 16.09. 2013. ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
// Variant 1: qS_u is included in the iteration factor in directions one and two, from k-th iteration k = 0,1,....
	
int Point::compute_pqr_explicit_full_v1 (int g, int m){
	
	// And ultimately this is the FULLY explicit version:
	
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	double margin = 1e-4;
	
	int N = Globals::getN();
	rays[g][m].dt_u = new (nothrow) double [N];
	
	int n;
	
	// Now distance in 3D space:
	double distance = sqrt ((X - rays[g][m].Xu) * (X - rays[g][m].Xu) + (Y - rays[g][m].Yu) * (Y - rays[g][m].Yu)) / Globals::getcostheta(g);
	double dt;
	
	for (n=0; n<N; n++)
		rays[g][m].dt_u[n] = distance * Globals::getprofile(n);
		
	// Here we compute p,q,r with the assumption of parabolic shape of the source function:
		
	for (n=0; n<N; n++){
		
		dt = rays[g][m].dt_u[n];
		
		if (dt < e-20){
			p[g][m][n] = 0.0;
			q[g][m][n] = 0.0;
			r[g][m][n] = 0.0;
			exxp[g][m][n] = 1.0;
		}
		
		else if (dt < margin){
			
			p[g][m][n] = 2.0/3.0 * dt  - dt*dt/4.0 + dt*dt*dt/15.0;
			q[g][m][n] = dt / 3.0 - dt * dt / 4.0 + dt * dt * dt / 10.0;
			r[g][m][n] = -dt * dt / 6.0 + dt * dt * dt / 12.0 - dt * dt * dt * dt/40;
			exxp[g][m][n] = 1.0 - dt + dt * dt * 0.5 - dt * dt * dt / 6.0;
		}
		
		else {
			
			p[g][m][n] = 1.0 + 2.0/dt * exp(-dt) - 2.0/dt/dt * (1.0 - exp(-dt));
			q[g][m][n] = -exp(-dt) - 2.0/dt * exp(-dt) + 2.0/dt/dt * (1.0 - exp(-dt));
			r[g][m][n] = -1.0 - exp(-dt) + 2.0/dt * (1.0 - exp(-dt));
			exxp[g][m][n] = exp(-dt);
		}
	}
// ------------------------------------------------------------------------------------------------------------------------------

	// But now some modifications have to be made:
	
	// Then	 we modify for the contribution of the local source function to the x derivative:
	
	for (n=0; n<N; n++){
			
		if (i==0){
		
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / (Globals::getX(0) - Globals::getX(1));
			
		
		}
		 
		else if (i == Globals::getNX() -1){
	
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / (Globals::getX(i) - Globals::getX(i-1));
			
		}
			
		else{
			
			double x1 = Globals::getX(i-1);
			double x2 = Globals::getX(i);
			double x3 = Globals::getX(i+1);
		
			p[g][m][n] += r[g][m][n] * Globals::getcosphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n)
				* 1.0 / (x3 - x1) * ((x3 - x2) / (x2 - x1) + (x1 - x2) / (x3 - x2));
			
		}
	}
	
		
// ---------------------------------------------------------------------------------------------------------------------------------
		
	// Then we modify for the contribution of the local source function to the y derivative:
	
	
	for (n=0; n<N; n++){
			
		if (j == 0){
		
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / ( Globals::getY(j) - Globals::getY(j+1) );
		
		}
			
		else if (j == Globals::getNY()-1){
		
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) / ( Globals::getY(j) - Globals::getY(j-1) );
		
		}
			
		else {
			
			double y1 = Globals::getY(j-1);
			double y2 = Globals::getY(j);
			double y3 = Globals::getY(j+1);
			
			p[g][m][n] += r[g][m][n] * Globals::getsinphi_C(g, m) * Globals::getcostheta(g) / Globals::getprofile(n) *
				1.0 / (y3 - y1) * ((y3 - y2) / (y2 - y1) + (y1 - y2) / (y3 - y2));
				
		}	
	}
	
	// And compute r1 .. r4 :
	
	double x1, x2, x3;
	double y1, y2, y3;
	double x, y;
		
	for (n=0; n<N; n++){
		
		// First r1 and r2, which arise due to the derivative with respect to x:
		
		if (i == 0){
			x1 = Globals::getX(i);
			x2 = Globals::getX(i+1);
			r2[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) / (x2 - x1);
		}
			
		else if (i == NX-1){
			x1 = Globals::getX(i-1);
			x2 = Globals::getX(i);
			r1[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) / -(x2 - x1);
		}
			
		else {
			x1 = Globals::getX(i-1);
			x2 = Globals::getX(i);
			x3 = Globals::getX(i+1);
			r1[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) *
			(x2 - x3) / (x1 - x2) / (x1 - x3);
			r2[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n) *
			(x2 - x1) / (x3 - x1) / (x3 - x2);
		} 
		
		// And then r3 and r4, which arise due to the derivative with respect to y:
		
		if (j == 0){
			y1 = Globals::getY(j);
			y2 = Globals::getY(j+1);
			r4[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) / (y2 - y1);
		}
			
		else if (j == NX-1){
			y1 = Globals::getY(j);
			y2 = Globals::getY(j+1);
			r3[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) / -(y2 - y1);
		}
		
		else {
			y1 = Globals::getY(j-1);
			y2 = Globals::getY(j);
			y3 = Globals::getY(j+1);
			r3[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) *
			(y2 - y3) / (y1 - y2) / (y1 - y3);
			r4[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n) *
			(y2 - y1) / (y3 - y1) / (y3 - y2);
		}	
	}
	
	// ----------------------------------------------------------------------------------------------------------------------------
	
	// And then modify them further due to the contribution of surrounding source functions to the upwind source function:
	
	// For turning on and off factors and contributions:
	
	double factor1, factor2, factor3, factor4;
	
	// 1: there is factor, 0 there is no factor.
	
	factor1 = 0.0;
	factor2 = 0.0;
	factor3 = 0.0;
	factor4 = 0.0;
	
	
	// QUADRANT 1:
	
	
	if (Globals::getphi_C(g, m) <= pi/2) {
		if (rays[g][m].type_u == 1){
			
			if (i == 0){
				for (n=0; n<N; n++){
					r5[g][m][n] += 0.0;
					r3[g][m][n] += 0.0;
					r6[g][m][n] += 0.0;
				}
			}
			
			else if (i == NX-1){
				
				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x = rays[g][m].Xu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += (1.0 - factor1) * (x2 - x) / (x2 - x1) * q[g][m][n];
					r3[g][m][n] += (1.0 - factor1) * (x1 - x) / (x1 - x2) * q[g][m][n];
					r6[g][m][n] += 0.0;
				}
			}
			
			else {
				
				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(i+1);
				x = rays[g][m].Xu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += (1.0 - factor1) * (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
					r3[g][m][n] += (1.0 - factor1) * (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
					r6[g][m][n] += (1.0 - factor1) * (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
				}
			}
			
		}
		
		else if (rays[g][m].type_u == 2){
			
			if (j == 0){
				for (n=0; n<N; n++){
					r5[g][m][n] += 0.0;
					r1[g][m][n] += 0.0;
					r7[g][m][n] += 0.0;
				}
			}
			
			else if (j == NY-1){
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += (1.0 - factor1) * (y2 - y) / (y2 - y1) * q[g][m][n];
					r1[g][m][n] += (1.0 - factor1) * (y1 - y) / (y1 - y2) * q[g][m][n];
					r7[g][m][n] += 0.0;
				}
			}
			
			else {
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y3 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += (1.0 - factor1) * (y2 - y) * (y3 - y) / (y2 - y1) / (y3 - y1) * q[g][m][n];
					r1[g][m][n] += (1.0 - factor1) * (y1 - y) * (y3 - y) / (y1 - y2) / (y3 - y2) * q[g][m][n];
					r7[g][m][n] += (1.0 - factor1) * (y1 - y) * (y2 - y) / (y1 - y3) / (y2 - y3) * q[g][m][n];
				}
			}
		}
	}
	
	// QUADRANT 2:
	
	else if (Globals::getphi_C(g, m) <= pi){
		if (rays[g][m].type_u == 1){
			
			if (i == NX - 1){
				for (n=0; n<N; n++){
					r5[g][m][n] += 0.0;
					r3[g][m][n] += 0.0;
					r6[g][m][n] += 0.0;
				}
			}
			
			else if (i == 0){
				
				x1 = Globals::getX(i);
				x2 = Globals::getX(i+1);
				x = rays[g][m].Xu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += 0.0;
					r3[g][m][n] += (1.0 - factor2) * (x2 - x) / (x2 - x1) * q[g][m][n];
					r6[g][m][n] += (1.0 - factor2) * (x1 - x) / (x1 - x2) * q[g][m][n];
				}
			}
			
			else {
				
				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(i+1);
				x = rays[g][m].Xu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += (1.0 - factor2) * (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
					r3[g][m][n] += (1.0 - factor2) * (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
					r6[g][m][n] += (1.0 - factor2) * (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
				}
			}
		}
		
		else if (rays[g][m].type_u == 2){
			
			if (j == 0){
				for (n=0; n<N; n++){
					r6[g][m][n] += 0.0;
					r2[g][m][n] += 0.0;
					r8[g][m][n] += 0.0;
				}
			}
			
			else if (j == NY-1){
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r6[g][m][n] += (1.0 - factor2) * (y2 - y) / (y2 - y1) * q[g][m][n];
					r2[g][m][n] += (1.0 - factor2) * (y1 - y) / (y1 - y2) * q[g][m][n];
					r8[g][m][n] += 0.0;
				}
			}
			
			else {
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y3 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r6[g][m][n] += (1.0 - factor2) * (y2 - y) * (y3 - y) / (y2 - y1) / (y3 - y1) * q[g][m][n];
					r2[g][m][n] += (1.0 - factor2) * (y1 - y) * (y3 - y) / (y1 - y2) / (y3 - y2) * q[g][m][n];
					r8[g][m][n] += (1.0 - factor2) * (y1 - y) * (y2 - y) / (y1 - y3) / (y2 - y3) * q[g][m][n];
				}
			}
		}
	}
	
	// QUADRANT 3:
	
	else if (Globals::getphi_C(g, m) <= pi * 3 / 2){
		if (rays[g][m].type_u == 1){
			
			if (i == NX - 1){
				for (n=0; n<N; n++){
					r7[g][m][n] += 0.0;
					r4[g][m][n] += 0.0;
					r8[g][m][n] += 0.0;
				}
			}
			
			else if (i == 0){
				
				x1 = Globals::getX(i);
				x2 = Globals::getX(i+1);
				x = rays[g][m].Xu;
				
				for (n=0; n<N; n++){
					r7[g][m][n] += 0.0;
					r4[g][m][n] += (x2 - x) / (x2 - x1) * q[g][m][n];
					r8[g][m][n] += (x1 - x) / (x1 - x2) * q[g][m][n];
				}
			}
			
			else {
				
				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(i+1);
				x = rays[g][m].Xu;
				
				for (n=0; n<N; n++){
					r7[g][m][n] += (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
					r4[g][m][n] += (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
					r8[g][m][n] += (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
				}
			}
		}
		
		else if (rays[g][m].type_u == 2){
			
			if (j == NY-1){
				for (n=0; n<N; n++){
					r6[g][m][n] += 0.0;
					r2[g][m][n] += 0.0;
					r8[g][m][n] += 0.0;
				}
			}
			
			else if (j == 0){
				
				y1 = Globals::getY(j);
				y2 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r6[g][m][n] += 0.0;
					r2[g][m][n] += (y2 - y) / (y2 - y1) * q[g][m][n];
					r8[g][m][n] += (y1 - y) / (y1 - y2) * q[g][m][n];
				}
				
			}
			
			else {
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y3 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r6[g][m][n] += (y2 - y) * (y3 - y) / (y2 - y1) / (y3 - y1) * q[g][m][n];
					r2[g][m][n] += (y1 - y) * (y3 - y) / (y1 - y2) / (y3 - y2) * q[g][m][n];
					r8[g][m][n] += (y1 - y) * (y2 - y) / (y1 - y3) / (y2 - y3) * q[g][m][n];
				}
			}
		}
	}
	
	// QUADRANT 4:
	
	else if (Globals::getphi_C(g, m) <= pi * 4 / 2){
		if (rays[g][m].type_u == 1){
			
			if (i == 0){
				for (n=0; n<N; n++){
					r7[g][m][n] += 0.0;
					r4[g][m][n] += 0.0;
					r8[g][m][n] += 0.0;
				}
			}
			
			else if (i == NX-1){
				
				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x = rays[g][m].Xu;
				
				for (n=0; n<N; n++){
					r7[g][m][n] += (x2 - x) / (x2 - x1) * q[g][m][n];
					r4[g][m][n] += (x1 - x) / (x1 - x2) * q[g][m][n];
					r8[g][m][n] += 0.0;
				}
			}
			
			else {
				
				x1 = Globals::getX(i-1);
				x2 = Globals::getX(i);
				x3 = Globals::getX(i+1);
				x = rays[g][m].Xu;
				
				for (n=0; n<N; n++){
					r7[g][m][n] += (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1) * q[g][m][n];
					r4[g][m][n] += (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2) * q[g][m][n];
					r8[g][m][n] += (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3) * q[g][m][n];
				}
			}
		}
		
		else if (rays[g][m].type_u == 2){
			
			if (j == NY-1){
				for (n=0; n<N; n++){
					r5[g][m][n] += 0.0;
					r1[g][m][n] += 0.0;
					r7[g][m][n] += 0.0;
				}
			}
			
			else if (j == 0){
				
				y1 = Globals::getY(j);
				y2 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += 0.0;
					r1[g][m][n] += (y2 - y) / (y2 - y1) * q[g][m][n];
					r7[g][m][n] += (y1 - y) / (y1 - y2) * q[g][m][n];
				}
			}
			
			else {
				
				y1 = Globals::getY(j-1);
				y2 = Globals::getY(j);
				y3 = Globals::getY(j+1);
				y = rays[g][m].Yu;
				
				for (n=0; n<N; n++){
					r5[g][m][n] += (y2 - y) * (y3 - y) / (y2 - y1) / (y3 - y1) * q[g][m][n];
					r1[g][m][n] += (y1 - y) * (y3 - y) / (y1 - y2) / (y3 - y2) * q[g][m][n];
					r7[g][m][n] += (y1 - y) * (y2 - y) / (y1 - y3) / (y2 - y3) * q[g][m][n];
				}
			}
		}
	}
	
	delete []rays[g][m].dt_u;
	
	return 0;
	
	
}

int Point::formal_solution_explicit_full_v1 (int g, int m, Point ** Grid, int iter){
	
	// First establish four neighboring source functions:
	
	double S1, S2, S3, S4, S5, S6, S7, S8;
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	// First "x" ones:
	
	if (i==0){
		S1 = 0;
		S2 = Grid[i+1][j].S;
	}
	
	else if (i== Globals::getNX() - 1){
		S1 = Grid[i-1][j].S;
		S2 = 0;
	}
	
	else {
		S1 = Grid[i-1][j].S;
		S2 = Grid[i+1][j].S;
	}
	
	
	// Now "y" ones:
	
	if (j==0){
		S3 = 0;
		S4 = Grid[i][j+1].S;
	}
	
	else if (j == Globals::getNY() - 1){
		S3 = Grid[i][j-1].S;
		S4 = 0;
	}
	
	else {
		S3 = Grid[i][j-1].S;
		S4 = Grid[i][j+1].S;
	}
	
	// And corner ones:
	
	if (i == 0 || j == 0)
		S5 = 0.0;
	
	else 
		S5 = Grid[i-1][j-1].S;

	if (i== 0 || j == NY-1)
		S7 = 0.0;
	else 
		S7 = Grid[i-1][j+1].S;
		
	if (i == NX - 1 || j == 0)
		S6 = 0.0;
	else 
		S6 = Grid[i+1][j-1].S;
		
	if (i == NX - 1 || j == NY-1)
		S8 = 0.0;
	else 
		S8 = Grid[i+1][j+1].S;
		
	// Compute intensity and add appropriate coefficients:
	
	double factor1 = 0.0;
	double factor2 = 0.0;
	double factor3 = 0.0;
	double factor4 = 0.0;
	
	for (int n = 0; n<Globals::getN(); n++){
		
		I[g][m][n] = rays[g][m].Iu[n] * exxp[g][m][n] + p[g][m][n] * S  + 
			r1[g][m][n] * S1 + r2[g][m][n] * S2 + r3[g][m][n] * S3 + r4[g][m][n] * S4 + 
			r5[g][m][n] * S5 + r6[g][m][n] * S6 + r7[g][m][n] * S7 + r8[g][m][n] * S8;
			
		if (m < Globals::getNP_C(g) * 1 / 4){
			AA1 += 0.0 * (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB1 += (p[g][m][n] + rays[g][m].Iu[n] * exxp[g][m][n] / S + factor1 * q[g][m][n] * rays[g][m].Su / S) *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 2 / 4){
			AA2 += 0.0 * (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB2 += (p[g][m][n] + rays[g][m].Iu[n] * exxp[g][m][n] / S + factor2 * q[g][m][n] * rays[g][m].Su / S) *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 3 / 4){
			AA3 += (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB3 += p[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 4 / 4){
			AA4 += (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB4 += p[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
			 
		if (iter == 0){
		
		C1 += r1[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C2 += r2[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C3 += r3[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C4 += r4[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C5 += r5[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C6 += r6[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C7 += r7[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		C8 += r8[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
	}
	
	delete []rays[g][m].Iu;
	
	return 0;
	
}

// Variant 2: Factors in three directions:
	
int Point::formal_solution_explicit_full_v2 (int g, int m, Point ** Grid, int iter){
	
	// First establish four neighboring source functions:
	
	double S1, S2, S3, S4, S5, S6, S7, S8;
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	// First "x" ones:
	
	if (i==0){
		S1 = 0;
		S2 = Grid[i+1][j].S;
	}
	
	else if (i== Globals::getNX() - 1){
		S1 = Grid[i-1][j].S;
		S2 = 0;
	}
	
	else {
		S1 = Grid[i-1][j].S;
		S2 = Grid[i+1][j].S;
	}
	
	
	// Now "y" ones:
	
	if (j==0){
		S3 = 0;
		S4 = Grid[i][j+1].S;
	}
	
	else if (j == Globals::getNY() - 1){
		S3 = Grid[i][j-1].S;
		S4 = 0;
	}
	
	else {
		S3 = Grid[i][j-1].S;
		S4 = Grid[i][j+1].S;
	}
	
	// And corner ones:
	
	if (i == 0 || j == 0)
		S5 = 0.0;
	
	else 
		S5 = Grid[i-1][j-1].S;

	if (i== 0 || j == NY-1)
		S7 = 0.0;
	else 
		S7 = Grid[i-1][j+1].S;
		
	if (i == NX - 1 || j == 0)
		S6 = 0.0;
	else 
		S6 = Grid[i+1][j-1].S;
		
	if (i == NX - 1 || j == NY-1)
		S8 = 0.0;
	else 
		S8 = Grid[i+1][j+1].S;
		
	// Compute intensity and add appropriate coefficients:
	
	double factor1 = 0.0;
	double factor2 = 0.0;
	double factor3 = 0.0;
	double factor4 = 0.0;
	
	for (int n = 0; n<Globals::getN(); n++){
		
		I[g][m][n] = rays[g][m].Iu[n] * exxp[g][m][n] + p[g][m][n] * S  + 
			r1[g][m][n] * S1 + r2[g][m][n] * S2 + r3[g][m][n] * S3 + r4[g][m][n] * S4 + 
			r5[g][m][n] * S5 + r6[g][m][n] * S6 + r7[g][m][n] * S7 + r8[g][m][n] * S8;
			
		if (m < Globals::getNP_C(g) * 1 / 4){
			AA1 += 0.0 * (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB1 += (p[g][m][n] + rays[g][m].Iu[n] * exxp[g][m][n] / S) *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 2 / 4){
			AA2 += 0.0 * (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB2 += (p[g][m][n] + rays[g][m].Iu[n] * exxp[g][m][n] / S) *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 3 / 4){
			AA3 += (1- (iter - iter/2*2)) * (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB3 += (p[g][m][n] +  (iter - iter/2*2) * rays[g][m].Iu[n] * exxp[g][m][n]/S) *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
		
		else if (m < Globals::getNP_C(g) * 4 / 4){
			AA4 += (iter - iter/2*2) * (rays[g][m].Iu[n] * exxp[g][m][n]) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			BB4 += (p[g][m][n] + (1 - iter + iter/2*2) * rays[g][m].Iu[n] * exxp[g][m][n]/S) *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		}
	}
	
	delete []rays[g][m].Iu;
	
	return 0;
	
}

// Another take on the classic variant

int Point::formal_solution_classic (int g, int m, Point ** Grid, int iter){
	
	// Here we compute the intensity:
	
	for (int n = 0; n<Globals::getN(); n++)
		I[g][m][n] = rays[g][m].Iu[n] * exxp[g][m][n] + p[g][m][n] * S + q[g][m][n] * rays[g][m].Su + rx[g][m][n] * Spx + ry[g][m][n] * Spy;
			
	// But in the heart of FBILI method is the computation of the AA and BB coefficients:
	
	// For the moment let's compute them Jacobi style:
	
	for (int n = 0; n<Globals::getN(); n++){
		
		AA += (rays[g][m].Iu[n] * exxp[g][m][n] + q[g][m][n] * rays[g][m].Su) * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		BB += p[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		CX += rx[g][m][n] * 
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		CY += ry[g][m][n] *
			0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			
		// But then BB has to be additionally modified:
		
		// After that we can compute the modifications to the local source function:
		
		/*if (i == 0)
			BB += rx[g][m][n] *  - 1.0 / (Globals::getX(i+1) - Globals::getX(i)) * 0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		
		else if (i  == 1)
			BB += rx[g][m][n] * 1.0 / (Globals::getX(i) - Globals::getX(i-1)) * 0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			
		else
			BB += rx[g][m][n] * 2.0 / (Globals::getX(i) - Globals::getX(i-1)) * 0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			
		
		
		if (j == Globals::getNY()-1)
			BB += ry[g][m][n] * +1.0 / (Globals::getY(j) - Globals::getY(j-1)) * 0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
		
		else if (j  == Globals::getNY()-2)
			BB += ry[g][m][n] * -1.0 / (Globals::getY(j+1) - Globals::getY(j)) * 0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);
			
		else
			BB += ry[g][m][n] * -2.0 / (Globals::getY(j+1) - Globals::getY(j)) * 0.25 * Globals::getprofile(n) * Globals::getwx(n) * Globals::getw_C(g,m);*/
	
			
	}
	
	// Delete the upwind intensity as it is no longer needed:
	
	delete []rays[g][m].Iu;
	
	return 0;
}


int Point::compute_pqr_classic (int g, int m){
	
	// This is a function which computes p,q,r coefficients with the assumption that Spy and Spz are symmetrically computed 
	// according to the last direction of sweeping the grid
	
	double margin = 1e-4;
	
	int N = Globals::getN();
	rays[g][m].dt_u = new (nothrow) double [N];
	
	int n;
	
	// Now distance in 3D space:
	double distance = sqrt ((X - rays[g][m].Xu) * (X - rays[g][m].Xu) + (Y - rays[g][m].Yu) * (Y - rays[g][m].Yu)) / Globals::getcostheta(g);
	double dt;
	
	for (n=0; n<N; n++)
		rays[g][m].dt_u[n] = distance * Globals::getprofile(n);
		
	// Here we compute p,q,r with the assumption of parabolic shape of the source function:
		
	for (n=0; n<N; n++){
		
		dt = rays[g][m].dt_u[n];
		
		if (dt < e-20){
			p[g][m][n] = 0.0;
			q[g][m][n] = 0.0;
			r[g][m][n] = 0.0;
			exxp[g][m][n] = 1.0;
		}
		
		else if (dt < margin){
			
			p[g][m][n] = 2.0/3.0 * dt  - dt*dt/4.0 + dt*dt*dt/15.0;
			q[g][m][n] = dt / 3.0 - dt * dt / 4.0 + dt * dt * dt / 10.0;
			r[g][m][n] = -dt * dt / 6.0 + dt * dt * dt / 12.0 - dt * dt * dt * dt/40;
			exxp[g][m][n] = 1.0 - dt + dt * dt * 0.5 - dt * dt * dt / 6.0;
		}
		
		else {
			
			p[g][m][n] = 1.0 + 2.0/dt * exp(-dt) - 2.0/dt/dt * (1.0 - exp(-dt));
			q[g][m][n] = -exp(-dt) - 2.0/dt * exp(-dt) + 2.0/dt/dt * (1.0 - exp(-dt));
			r[g][m][n] = -1.0 - exp(-dt) + 2.0/dt * (1.0 - exp(-dt));
			exxp[g][m][n] = exp(-dt);
		}
	}
// ------------------------------------------------------------------------------------------------------------------------------

	// And compute rx, and ry:
		
	for (n=0; n<N; n++){
			
		rx[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getcosphi_C(g, m) / Globals::getprofile(n);
		ry[g][m][n] = r[g][m][n] * Globals::getcostheta(g) * Globals::getsinphi_C(g, m) / Globals::getprofile(n);
			
	}
	
	
	delete []rays[g][m].dt_u;
	
	return 0;
}

int Point::compute_S_classic (Point ** Grid){
	
	double Sp_x_mod, Sp_y_mod;
	
	// First modify BB because of Cx and Cy, and compute so-called modified derivatives which go into the expression for the source function
	/*
	if (i == 0 && j == Globals::getNY()-1){
		BB += 0.0;
		Sp_x_mod = 0.0;
		Sp_y_mod = 0.0;
	}
	
	else if (i == 0){
		BB += - CY / (Globals::getY(j+1) - Globals::getY(j));
		Sp_x_mod = 0.0;
		Sp_y_mod = Grid[i][j+1].S / (Globals::getY(j+1) - Globals::getY(j));
	}
	
	else if (j == Globals::getNY() - 1){
		BB += CX / (Globals::getX(i) - Globals::getx(i-1));
		Sp_x_mod = - Grid[i-1][i].S / (Globals::getX(i) - Globals::getx(i-1));
		Sp_y_mod = 0.0;
	}
		
	else if (j == Globals::getNY() - 2 && i>=1){
		BB += - CY / (Globals::getY(j+1) - Globals::getY(j)) + CX / (Globals::getX(i) - Globals::getx(i-1));
		Sp_x_mod = - Grid[i-1][i].S / (Globals::getX(i) - Globals::getx(i-1));
		Sp_y_mod = Grid[i][j+1].S / (Globals::getY(j+1) - Globals::getY(j));
	}
	
	else if (i == 1 && j < Globals::getNY()-1){
		BB += - CY / (Globals::getY(j+1) - Globals::getY(j)) + CX / (Globals::getX(i) - Globals::getx(i-1));
		Sp_x_mod = - Grid[i-1][i].S / (Globals::getX(i) - Globals::getx(i-1));
		Sp_y_mod = Grid[i][j+1].S / (Globals::getY(j+1) - Globals::getY(j));
	}
		
	else {
		BB +=  - 2.0 * CY / (Globals::getY(j+1) - Globals::getY(j)) + 2.0 * CX / (Globals::getX(i) - Globals::getx(i-1));
		Sp_x_mod = - 2.0 * Grid[i-1][i].S / (Globals::getX(i) - Globals::getx(i-1)) - Grid[i-1][j].Spx;
		Sp_y_mod = 2.0 * Grid[i][j+1].S / (Globals::getY(j+1) - Globals::getY(j)) - Grid[i][j+1].Spy;
	}
	*/
	// Alternate variant (how I think it should be done):
	
	if (i ==0 ){
		BB += - CX / (Globals::getX(i+1) - Globals::getX(i));
		Sp_x_mod = Grid[i+1][j].S / (Globals::getX(i+1) - Globals::getX(i));
	}
	else if (i == 1){
		BB += CX / (Globals::getX(i) - Globals::getX(i-1));
		Sp_x_mod = -Grid[i-1][j].S / (Globals::getX(i) - Globals::getX(i-1));
	}
	else {
		BB += 2.0 * CX / (Globals::getX(i) - Globals::getX(i-1));
		Sp_x_mod = -2.0 * Grid[i-1][j].S / (Globals::getX(i) - Globals::getX(i-1)) - Grid[i-1][j].Spx;
	}
	
	if (j == Globals::getNY()-1 ){
		BB += CY / (Globals::getY(j) - Globals::getY(j-1));
		Sp_y_mod = -  Grid[i][j-1].S / (Globals::getY(j) - Globals::getY(j-1));
	}
	else if (j == Globals::getNY()-2){
		BB += -CY / (Globals::getY(j+1) - Globals::getY(j));
		Sp_y_mod = Grid[i][j+1].S / (Globals::getY(j+1) - Globals::getY(j));
	}
	else {
		BB += - 2.0 * CY / (Globals::getY(j+1) - Globals::getY(j));
		Sp_y_mod = 2.0 * Grid[i][j+1].S / (Globals::getY(j+1) - Globals::getY(j)) - Grid[i][j+1].Spy;
	}
		
	
	// Update coefficient A accordingly
		
	AA += CX * Sp_x_mod + CY * Sp_y_mod;
	
	// Compute change in the source function and update the source function
	
	dS = (((1-ep) * AA + ep * B)) / (1 - (1-ep) * BB) - S;
	 
	S += dS;
	
	// Update the derivative
	
	if (i==0)
		Spx = 0;
	else 
		Spx = (Grid[i][j].S - Grid[i-1][j].S) / (Globals::getX(i) - Globals::getX(i-1));
		
	if (j==Globals::getNY()-1)
		Spy = 0;
	else 
		Spy = (Grid[i][j+1].S - Grid[i][j].S) / (Globals::getY(j+1) - Globals::getY(j));
	
	for (int g=0; g<Globals::getNT(); g++)
		for (int m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++)
			for (int n=0; n<Globals::getN(); n++)
				rays[g][m].Iu[n] * exxp[g][m][n] + p[g][m][n] * S + q[g][m][n] * rays[g][m].Su + rx[g][m][n] * Spx + ry[g][m][n] * Spy;
	
	return 0;
}

// --------------------------------------------------------------------------------------------------------------------------------------


