#include "all.h"

using namespace std;

int Globals::N = 0;
double * Globals::x;// = new (nothrow) double [0];
double * Globals::wx;// = new (nothrow) double [0];
double * Globals::profile;// = new (nothrow) double [0];
int Globals::type = 1;

int Globals::NX = 0;
int Globals::NY = 0;

int Globals::NP = 0;
int Globals::NT = 0;
int * Globals::NP_Carlson;// = new (nothrow) int [0];
double * Globals::phi;// = new (nothrow) double [0];
double ** Globals::phi_C;// = new (nothrow) double * [0];
double ** Globals::w_C;// = new (nothrow) double * [0];
double * Globals::theta;// = new (nothrow) double [0];
double * Globals::sintheta;// = new (nothrow) double [0];
double * Globals::costheta;// = new (nothrow) double [0];
double * Globals::X;// = new (nothrow) double [0];
double * Globals::Y;// = new (nothrow) double [0];

double ** Globals::sinphi_C;// = new (nothrow) double * [0];
double ** Globals::cosphi_C;// = new (nothrow) double * [0];

double **** Globals::Psi;

double Globals::w = 1;
int Globals::factor = 0;

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Functions:

// Profile-related:

void Globals::allocate_x (int N_in){
	N = N_in;
	x = new (nothrow) double [N];
}

void Globals::compute_x (){
	int n;

	if (N == 1)
		x[0] = 0;

	else if (type == 1){ // Doppler grid
		for (n=0; n<N; n++)
			x[n] = n * 5.0/(N-1);
	}

	else if (type == 2){ // Voigt grid
		double core = 4;
		int N_core = core / 0.5 + 1;
		int N_wing = N - N_core;
		for (n=0; n<N_core; n++)
			x[n] = n * 0.5;
		for (n=N_core; n<N; n++)
			x[n] = x[n-1] + pow(10, 0.1 * (n - N_core)); 
	}
}

double Globals::getx(int n){
	return x[n];
}
	
void Globals::compute_profile(){

	int n;
	profile = new (nothrow) double [N];
	wx = new (nothrow) double [N];

	if (N == 1){
		profile[0] = 1.0;
		wx[0] = 1.0;
	}

	else{ 

		for (n=0; n<N; n++){
			if (type == 1) profile[n] = 1.0/pisqrt * pow (e, - x[n] * x[n]);
			else if (type == 2) profile[n] = voigt(0.001, x[n]);
		}

		wx[0] = 0.5 * (x[1] - x[0]);
		wx[N-1] = 0.5 * (x[N-1] - x[N-2]);
		for (n=1; n<N-1; n++)
			wx[n] = 0.5 * (x[n+1] - x[n-1]);

		// Now a small chapter to renormalize weights if needed:
		double norm = 0;
		for (n=0; n<N; n++)
			norm += profile[n] * wx[n];
		for (n=0; n<N; n++)
			wx[n] /= norm;
	}
}
	
double Globals::getprofile(int n){
	return profile[n];
}

double Globals::getwx(int n){
	return wx[n];
}

int Globals::getN(){
	return N;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void Globals::setNX(int NXin){
	NX = NXin;
}

int Globals::getNX(){
	return NX;
}

void Globals::setNY(int NYin){
	NY = NYin;
}

int Globals::getNY(){
	return NY;
}


// ---------------------------------------------------------------------------
// Now a section for angles for cylindrical geometry:

void Globals::setNP (int NP_in){
	NP = NP_in;
}
	
void Globals::setNT (int NT_in){
	NT = NT_in;
}

int Globals::getNP (){
	return NP;
}
	
int Globals::getNT (){
	return NT;
}

double Globals::gettheta(int g){
	return theta[g];
}

double Globals::getsintheta (int g){
	return sintheta[g];
}

double Globals::getcostheta (int g){
	return costheta[g];
}

// And, special mention here, additional function which will deal with the angles according to article by Carlson (1963)

int Globals::compute_angles_Carlson (){
	
	// First, this grid is set-up (at the moment) with pretty much fixed number of angles.
	
	NT = 8;
	delete []NP_Carlson;
	NP_Carlson = new (nothrow) int [NT];
	NP_Carlson[0] = 4;
	NP_Carlson[1] = 6;
	NP_Carlson[2] = 8;
	NP_Carlson[3] = 6;
	NP_Carlson[4] = 6;
	NP_Carlson[5] = 8;
	NP_Carlson[6] = 6;
	NP_Carlson[7] = 4;
	
	double mu_l[4];
	mu_l[0] = 0.11104445;
	mu_l[1] = 0.50307327;
	mu_l[2] = 0.70273364;
	mu_l[3] = 0.85708017;
	
	int g, m;
	
	// Add the region from pi to 2pi
	
	for (g=0; g<NT; g++)
		NP_Carlson[g] *= 2;
	
	// Let us first compute theta:
	
	theta = new (nothrow) double [NT];
	
	for (g = 0; g<NT/2; g++)
		theta[4+g] = asin(mu_l[g]);
	
	for (g = 0; g<NT/2; g++)
		theta[g] = -theta[NT-1-g];
	
	// And now take care of the angle phi:
	
	phi_C = new (nothrow) double * [NT];
	
	for (g=0; g<NT; g++)
		* (phi_C + g) = new (nothrow) double [NP_Carlson[g]];
		
	// Now we by hand input all the values of mu_y or if you prefer mu_m in the grid. Only for the octant with positive theta and positive phi
	
	phi_C[4][0] = mu_l[3];
	phi_C[4][1] = mu_l[2];
	phi_C[4][2] = mu_l[1];
	
	phi_C[5][0] = mu_l[3];
	phi_C[5][1] = mu_l[2];
	phi_C[5][2] = mu_l[1];
	phi_C[5][3] = mu_l[0];
	
	phi_C[6][0] = mu_l[2];
	phi_C[6][1] = mu_l[1];
	phi_C[6][2] = mu_l[0];
	
	phi_C[7][0] = mu_l[1];
	phi_C[7][1] = mu_l[0];
	
	// Now we manualy modify them:
	
	for (g=NT/2; g<NT; g++)
		for (m = 0; m<NP_Carlson[g]/4; m++){
			phi_C[g][m] = acos (phi_C[g][m] / cos(theta[g]));
		}
		
	// And now it is time to assign integration weights:
	
	delete []w_C;
	w_C = new (nothrow) double * [NT];
	for (g=0; g<NT; g++)
		*(w_C + g) = new (nothrow) double [NP_Carlson[g]];
	
	w_C[4][0] = 0.11434821;
	w_C[4][2] = 0.11434821;
	w_C[5][0] = 0.11434821;
	w_C[5][3] = 0.11434821;
	w_C[7][0] = 0.11434821;
	w_C[7][1] = 0.11434821;
	
	w_C[4][1] = 0.07937696;
	w_C[6][0] = 0.07937696;
	w_C[6][2] = 0.07937696;
	
	w_C[5][1] = 0.02424996;
	w_C[5][2] = 0.02424996;
	w_C[6][1] = 0.02424996;
	
	// Now also for negative theta
		
	for (g=0; g<NT/2; g++)
		for (m=0; m<NP_Carlson[g]/4; m++){
			phi_C[g][m] = phi_C[NT-1-g][m];
			w_C[g][m] = w_C[NT-1-g][m];
		}
		
	// Now also for phi which is between pi/2 and pi:
	
	for (g=0; g<NT; g++)
		for (m=NP_Carlson[g]/4; m<NP_Carlson[g]/2; m++){
			phi_C[g][m] = pi - phi_C[g][NP_Carlson[g]/2- 1 - m];
			w_C[g][m] = w_C[g][NP_Carlson[g]/2 - 1 - m];
		}
		
	// And now all that for phi which is between pi and 2pi
	
	for (g=0; g<NT; g++)
		for (m=NP_Carlson[g]/2; m < NP_Carlson[g]; m++){
			phi_C[g][m] = 2.0 * pi - phi_C[g][NP_Carlson[g]-1-m];
			w_C[g][m] = w_C[g][NP_Carlson[g]-1-m];
			
		}
		
	// And test how it looks like:
	
	ofstream anglez;
	anglez.open("angles.txt");
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			anglez << g << "\t" << theta[g] << "\t" << m << "\t" << phi_C[g][m] << "\t" << w_C[g][m] << endl;
	anglez.close(); 
	
	// Now compute cosines and sines:
	
	delete []sinphi_C;
	delete []cosphi_C;
	
	sinphi_C = new (nothrow) double * [NT];
	cosphi_C = new (nothrow) double * [NT];
	
	for (g=0; g<NT; g++){
		* (sinphi_C + g) = new (nothrow) double [NP_Carlson[g]];
		* (cosphi_C + g) = new (nothrow) double [NP_Carlson[g]];
	}
	
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++){
			sinphi_C[g][m] = sin(phi_C[g][m]);
			cosphi_C[g][m] = cos(phi_C[g][m]);
		}
		
	// And also cosines and sines for theta:
	
	delete []sintheta;
	delete []costheta;
	
	sintheta = new (nothrow) double [NT];
	costheta = new (nothrow) double [NT];
	
	for (g=0; g<NT; g++){
		sintheta[g] = sin(theta[g]);
		costheta[g] = cos(theta[g]);
	}
	
	// Here we reduce theta to half range:
	
	NT /= 2;
	
	// To test the integration:
	
	double norma = 0;
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			norma += 0.25 * w_C[g][m];
	cout << "Norm is: " << norma << endl;
	
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			w_C[g][m] /= norma;
			
	return 0;
}

int Globals::compute_angles_Carlson_10 (){
	
	// We follow the same, manual, approach from the version with NT = 8.
	
	// First, this grid is set-up (at the moment) with pretty much fixed number of angles.
	
	int g, m;
	
	NT = 10;
	delete []NP_Carlson;
	NP_Carlson = new (nothrow) int [NT];
	NP_Carlson[0] = 4;
	NP_Carlson[1] = 6;
	NP_Carlson[2] = 8;
	NP_Carlson[3] = 10;
	NP_Carlson[4] = 8;
	NP_Carlson[5] = 8;
	NP_Carlson[6] = 10;
	NP_Carlson[7] = 8;
	NP_Carlson[8] = 6;
	NP_Carlson[9] = 4;
	
	for (g=0; g<NT; g++)
		NP_Carlson[g] *= 2;
	
	double mu_l[4];
	mu_l[0] = 0.10480667;
	mu_l[1] = 0.45209930;
	mu_l[2] = 0.63071635;
	mu_l[3] = 0.76890342;
	mu_l[4] = 0.88578880;
	
	// Let us first compute theta:
	
	theta = new (nothrow) double [NT];
	
	for (g = 0; g<NT/2; g++)
		theta[NT/2+g] = asin(mu_l[g]);
	
	for (g = 0; g<NT/2; g++)
		theta[g] = -theta[NT-1-g];	
	
	// And now take care of the angle phi:
	
	phi_C = new (nothrow) double * [NT];
	
	for (g=0; g<NT; g++)
		* (phi_C + g) = new (nothrow) double [NP_Carlson[g]];
		
	// Now we by hand input all the values of mu_y or if you prefer mu_m in the grid. Only for the octant with positive theta and positive phi
	
	phi_C[5][0] = mu_l[4];
	phi_C[5][1] = mu_l[3];
	phi_C[5][2] = mu_l[2];
	phi_C[5][3] = mu_l[1];
	
	phi_C[6][0] = mu_l[4];
	phi_C[6][1] = mu_l[3];
	phi_C[6][2] = mu_l[2];
	phi_C[6][3] = mu_l[1];
	phi_C[6][4] = mu_l[0];
	
	phi_C[7][0] = mu_l[3];
	phi_C[7][1] = mu_l[2];
	phi_C[7][2] = mu_l[1];
	phi_C[7][3] = mu_l[0];
	
	phi_C[8][0] = mu_l[2];
	phi_C[8][1] = mu_l[1];
	phi_C[8][2] = mu_l[0];
	
	phi_C[9][0] = mu_l[1];
	phi_C[9][1] = mu_l[0];
	
	// Now we manualy modify them:
	
	for (g=NT/2; g<NT; g++)
		for (m = 0; m<NP_Carlson[g]/4; m++){
			phi_C[g][m] = acos (phi_C[g][m] / cos(theta[g]));
		}
		
	// And now it is time to assign integration weights:
	
	delete []w_C;
	w_C = new (nothrow) double * [NT];
	for (g=0; g<NT; g++)
		*(w_C + g) = new (nothrow) double [NP_Carlson[g]];
	
	w_C[5][0] = 0.087396525;
	w_C[5][3] = 0.087396525;
	w_C[6][0] = 0.087396525;
	w_C[6][4] = 0.087396525;
	w_C[9][0] = 0.087396525;
	w_C[9][1] = 0.087396525;
	
	w_C[5][1] = 0.054876550;
	w_C[5][2] = 0.054876550;
	w_C[7][0] = 0.054876550;
	w_C[7][3] = 0.054876550;
	w_C[8][0] = 0.054876550;
	w_C[8][2] = 0.054876550;
	
	w_C[6][1] = 0.021936360;
	w_C[6][3] = 0.021936360;
	w_C[8][1] = 0.021936360;
	
	w_C[6][2] = 0.026850829;
	w_C[7][1] = 0.026850829;
	w_C[7][2] = 0.026850829;
	
	// Now also for negative theta
		
	for (g=0; g<NT/2; g++)
		for (m=0; m<NP_Carlson[g]/4; m++){
			phi_C[g][m] = phi_C[NT-1-g][m];
			w_C[g][m] = w_C[NT-1-g][m];
		}
		
	// Now also for phi which is between pi/2 and pi:
	
	for (g=0; g<NT; g++)
		for (m=NP_Carlson[g]/4; m<NP_Carlson[g]/2; m++){
			phi_C[g][m] = pi - phi_C[g][NP_Carlson[g]/2 -1 - m];
			w_C[g][m] = w_C[g][NP_Carlson[g]/2 -1 - m];
		}
	
	// And now, finally for the angles phi which fall between pi and 2pi
	
	for (g=0; g<NT; g++)
		for (m=NP_Carlson[g]/2; m<NP_Carlson[g]; m++){
			phi_C[g][m] = 2*pi - phi_C[g][NP_Carlson[g] -1 - m];
			w_C[g][m] = w_C[g][NP_Carlson[g] -1 - m];
		}
		
	// And test how it looks like:
	
	ofstream anglez;
	anglez.open("angles.txt");
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			anglez << g << "\t" << theta[g] << "\t" << m << "\t" << phi_C[g][m] << "\t" << w_C[g][m] << endl;
	anglez.close(); 
	
	// Now compute cosines and sines:
	
	delete []sinphi_C;
	delete []cosphi_C;
	
	sinphi_C = new (nothrow) double * [NT];
	cosphi_C = new (nothrow) double * [NT];
	
	for (g=0; g<NT; g++){
		* (sinphi_C + g) = new (nothrow) double [NP_Carlson[g]];
		* (cosphi_C + g) = new (nothrow) double [NP_Carlson[g]];
	}
	
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++){
			sinphi_C[g][m] = sin(phi_C[g][m]);
			cosphi_C[g][m] = cos(phi_C[g][m]);
		}
		
	// And also cosines and sines for theta:
	
	delete []sintheta;
	delete []costheta;
	
	sintheta = new (nothrow) double [NT];
	costheta = new (nothrow) double [NT];
	
	for (g=0; g<NT; g++){
		sintheta[g] = sin(theta[g]);
		costheta[g] = cos(theta[g]);
	}
	
	// To test the integration:
	
	double norma = 0;
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			norma += 0.25 * w_C[g][m];
	cout << "Norm is: " << norma << endl;
	
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			w_C[g][m] /= norma;
	

	return 0;
}

// And, another special mention here, we make a function which makes simplest possible quadrature : equidistant steps in theta and phi:

int Globals::compute_angles_simple (int NTnew, int NPnew){
	
	int g,m;
	
	// Free previous computations:
// ---------------------------------------------------------------------------------------------------------
	
	
	delete []theta;
	for (g=0; g<NT; g++){
		delete []phi_C[g];
		delete []w_C[g];
	}
		
	delete []phi_C;
	delete []w_C;
	
	delete []costheta;
	delete []sintheta;
	
	for (g=0; g<NT; g++){
		delete []cosphi_C[g];
		delete []sinphi_C[g];
	}
	
	delete []cosphi_C;
	delete []sinphi_C;
	
	delete []NP_Carlson; // Ignore this "Carlson" it is just the relict of the past in the name :)
// ---------------------------------------------------------------------------------------------------------
	
	NT = NTnew;
	NP_Carlson = new (nothrow) int [NT];

	for (g =0; g<NT; g++)
		NP_Carlson[g] = NPnew;
	
	theta = new (nothrow) double [NT];
	phi_C = new (nothrow) double * [NT]; 
	w_C = new (nothrow) double * [NT];
	
	
	for (g = 0; g<NT; g++){
		phi_C[g] = new (nothrow) double [NP_Carlson[g]];
		w_C[g] = new (nothrow) double [NP_Carlson[g]];
	}
	
	// We contruct simplest possible, equidistant grid:
	
	double steptheta = pi / (NT-1);
	double stepphi;
	
	for (g=0; g<NT; g++)
		theta[g] = - pi/2.0 + g * steptheta;
		
	theta[0] += 1e-3;
	theta[NT-1] -= 1e-3;
	
	for (g=0; g<NT; g++){
		
		stepphi = pi / (NP_Carlson[g] - 1);
		
		for (m=0; m<NP_Carlson[g]; m++)
			phi_C[g][m] = m * stepphi;
			
		
		phi_C[g][0] += 1e-3;
		phi_C[g][NP_Carlson[g]-1] -= 1e-3;
	}
	
	// Those were the angles, now the weights:
	
	double wt [NT];
	double wp [NP_Carlson[0]];
	
	for (g=0; g<NT; g++){
		if (g==0)
			wt[g] = theta[1] - theta[0];
		else if (g == NT-1)
			wt[g] = theta[NT-1] - theta[NT-2];	
		else 
			wt[g] = theta[g+1] - theta[g-1];
	}
	
	for (m=0; m<NP_Carlson[0]; m++){
		if (m == 0)
			wp[m] = phi_C[0][1] - phi_C[0][0];
		else if (m == NP_Carlson[0] -1)
			wp[m] = phi_C[0][m] - phi_C[0][m-1];
		else 
			wp[m] = phi_C[0][m+1] - phi_C[0][m-1];
	}
	
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			w_C[g][m] = cos(theta[g]) * wt[g] * wp[m];
			
	// Normalization:
// --------------------------------------------------------------------------------------------
	
	double norma = 0;
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			norma += 0.25 * w_C[g][m];
	cout << "Norm is: " << norma << endl;
	
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			w_C[g][m] /= norma;
			
	// Sines and cosines:
// ---------------------------------------------------------------------------------------------
	
	sintheta = new (nothrow) double [NT];
	costheta = new (nothrow) double [NT];
	sinphi_C = new (nothrow) double * [NT];
	cosphi_C = new (nothrow) double * [NT];
	
	for (g=0; g<NT; g++){
		sinphi_C[g] = new (nothrow) double [NP_Carlson[g]];
		cosphi_C[g] = new (nothrow) double [NP_Carlson[g]];
		
		sintheta[g] = sin(theta[g]);
		costheta[g] = cos(theta[g]);
		
		for (m=0; m<NP_Carlson[g]; m++){
			sinphi_C[g][m] = sin(phi_C[g][m]);
			cosphi_C[g][m] = cos(phi_C[g][m]);
		}
	}
	
// ---------------------------------------------------------------------------------------------
	cout << "New, simple, angular grid has been set-up" << endl;

	return 0;		
}
	
	
// And yet another special version of the compute angles, we call it "Pseudogaussian":

int Globals::compute_angles_pseudogauss (int NTnew, int NPnew){
	
	int g,m;
	
	// Free previous computations:
// ---------------------------------------------------------------------------------------------------------
	
	
	delete []theta;
	for (g=0; g<NT; g++){
		delete []phi_C[g];
		delete []w_C[g];
	}
		
	delete []phi_C;
	delete []w_C;
	
	delete []costheta;
	delete []sintheta;
	
	for (g=0; g<NT; g++){
		delete []cosphi_C[g];
		delete []sinphi_C[g];
	}
	
	delete []cosphi_C;
	delete []sinphi_C;
	
	delete []NP_Carlson; // Ignore this "Carlson" it is just the relict of the past in the name :)
// ---------------------------------------------------------------------------------------------------------
	
	// We allocate the memory first:
	
	NT = NTnew;
	NP_Carlson = new (nothrow) int [NT];

	for (g =0; g<NT; g++)
		NP_Carlson[g] = NPnew;
	
	theta = new (nothrow) double [NT];
	phi_C = new (nothrow) double * [NT]; 
	w_C = new (nothrow) double * [NT];
	
	
	for (g = 0; g<NT; g++){
		phi_C[g] = new (nothrow) double [NP_Carlson[g]];
		w_C[g] = new (nothrow) double [NP_Carlson[g]];
	}
	
// ----------------------------------------------------------------------------------------------------------
	
	double stepphi;
	double steptheta;
	
	alglib::real_1d_array mi;
	alglib::real_1d_array miw;
	alglib::ae_int_t info;
	
	mi.setlength(NT);
	miw.setlength(NT);
	
	alglib::gqgenerategausslegendre (NT, info, mi, miw);
	
	double * mitemp = mi.getcontent();
	double * miwtemp = miw.getcontent();
	
	for (g=0; g<NT; g++)
		theta[g] = asin(mitemp[g]);
		
	
	for (g=0; g<NT; g++){
		
		stepphi = 2.0 * pi / (NP_Carlson[g]);
		
		for (m=0; m<NP_Carlson[g]; m++)
			phi_C[g][m] = stepphi/2.0 + m * stepphi;
			
		
		phi_C[g][0] += 0;
		phi_C[g][NP_Carlson[g]-1] -= 0;
	}
	
	// Those were the angles, now the weights:
	
	double wt [NT];
	double wp [NP_Carlson[0]];
	
	for (g=0; g<NT; g++){
		if (g==0)
			wt[g] = theta[1] - theta[0];
		else if (g == NT-1)
			wt[g] = theta[NT-1] - theta[NT-2];	
		else 
			wt[g] = theta[g+1] - theta[g-1];
	}
	
	for (m=0; m<NP_Carlson[0]; m++){
		if (m == 0)
			wp[m] = phi_C[0][1] - phi_C[0][0];
		else if (m == NP_Carlson[0] -1)
			wp[m] = phi_C[0][m] - phi_C[0][m-1];
		else 
			wp[m] = phi_C[0][m+1] - phi_C[0][m-1];
			
		wp[m] = stepphi;
	}
	
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			w_C[g][m] = miwtemp[g] * wp[m];
					
	// Normalization:
// --------------------------------------------------------------------------------------------
	
	double norma = 0;
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			norma += 0.25 * w_C[g][m];
	cout << "Norm is: " << norma << endl;
	
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			w_C[g][m] /= norma;
			
	// Sines and cosines:
// ---------------------------------------------------------------------------------------------
	
	sintheta = new (nothrow) double [NT];
	costheta = new (nothrow) double [NT];
	sinphi_C = new (nothrow) double * [NT];
	cosphi_C = new (nothrow) double * [NT];
	
	for (g=0; g<NT; g++){
		sinphi_C[g] = new (nothrow) double [NP_Carlson[g]];
		cosphi_C[g] = new (nothrow) double [NP_Carlson[g]];
		
		sintheta[g] = sin(theta[g]);
		costheta[g] = cos(theta[g]);
		
		for (m=0; m<NP_Carlson[g]; m++){
			sinphi_C[g][m] = sin(phi_C[g][m]);
			cosphi_C[g][m] = cos(phi_C[g][m]);
		}
	}
	
// ---------------------------------------------------------------------------------------------
	cout << "New, simple, angular grid has been set-up" << endl;

	return 0;		
}


int Globals::compute_angles_1 (){
	
	// First, this grid is set-up (at the moment) with pretty much fixed number of angles.
	
	NT = 1;
	delete []NP_Carlson;
	NP_Carlson = new (nothrow) int [NT];
	NP_Carlson[0] = 4;
	
	int g, m;
	
	theta = new (nothrow) double [NT];
	
	theta[0] = - acos(sqrt(0.66666)); 
	
	// And now take care of the angle phi:
	
	phi_C = new (nothrow) double * [NT];
	w_C = new (nothrow) double * [NT];

	for (g=0; g<NT; g++){
		* (phi_C + g) = new (nothrow) double [NP_Carlson[g]];
		w_C[g] = new (nothrow) double [NP_Carlson[g]];
	}
	

	for (m=0; m<NP_Carlson[0]; m++){
		phi_C[0][m] = pi/2 + pi/2 * m;
		w_C[0][m] = 1.0;
	}

	//cout << "So far so good." << endl;	
	
		// And now it is time to assign integration weights:
		
	ofstream anglez;
	anglez.open("angles.txt");
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			anglez << g << "\t" << theta[g] << "\t" << m << "\t" << phi_C[g][m] << "\t" << w_C[g][m] << endl;
	anglez.close(); 
	
	// Now compute cosines and sines:
	
	delete []sinphi_C;
	delete []cosphi_C;
	
	sinphi_C = new (nothrow) double * [NT];
	cosphi_C = new (nothrow) double * [NT];
	
	for (g=0; g<NT; g++){
		* (sinphi_C + g) = new (nothrow) double [NP_Carlson[g]];
		* (cosphi_C + g) = new (nothrow) double [NP_Carlson[g]];
	}
	
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++){
			sinphi_C[g][m] = sin(phi_C[g][m]);
			cosphi_C[g][m] = cos(phi_C[g][m]);
		}
		
	// And also cosines and sines for theta:
	
	delete []sintheta;
	delete []costheta;
	
	sintheta = new (nothrow) double [NT];
	costheta = new (nothrow) double [NT];
	
	for (g=0; g<NT; g++){
		sintheta[g] = sin(theta[g]);
		costheta[g] = cos(theta[g]);
	}
		
	// To test the integration:
	
	double norma = 0;
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			norma += 0.25 * w_C[g][m];

	cout << "The simplest angular grid with one point per octant is used. Norm is: " << norma << endl;
	
	for (g=0; g<NT; g++)
		for (m=0; m<NP_Carlson[g]; m++)
			w_C[g][m] /= norma;
			
	return 0;
}


		
// Appropriate "get" functions:

double Globals::getphi_C (int g, int m){
	
	return phi_C[g][m];
}

double Globals::getcosphi_C (int g, int m){
	
	return cosphi_C[g][m];
}

double Globals::getsinphi_C (int g, int m){
		
	return sinphi_C[g][m];
}
	
double Globals::getw_C (int g, int m){
	
	return w_C[g][m];
}

int Globals::getNP_C (int g){
	
	return NP_Carlson[g];
}

// More geometrical stuff:

int Globals::set_geometry (Point ** Grid){

	int i, j;
	X = new (nothrow) double [NX];
	for (i=0; i<NX; i++)
		X[i] = Grid[i][0].X;

	Y = new (nothrow) double [NY];
	for (j=0; j<NY; j++)
		Y[j] = Grid[0][j].Y;

	return 0;
}
	
double * Globals::get_X_axis (){
	return X;
}

double * Globals::get_Y_axis (){
	return Y;
}

double Globals::getX (int i){
	return X[i];
}

double Globals::getY (int j){
	return Y[j];
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Deallocation of everything (DOES NOT WORK, SHOULD BE REPAIRED AT SOME POINT)

void Globals::deallocate_all(){
	
	delete []theta;
	delete []x;
	delete []profile;
	delete []wx;

	//delete []sinphi;
	//delete []cosphi;
	delete []sintheta;
	delete []costheta;
	
	int g;
	
	for (g=0; g<NT; g++){
		delete []phi_C[g];
		delete []sinphi_C[g];
		delete []cosphi_C[g];
		delete []w_C[g];
	}
	
	delete []phi_C;
	delete []sinphi_C;
	delete []cosphi_C;
	delete []w_C;
	
	delete []NP_Carlson;
	//for (m=0; m<NP; m++)
		//delete []jacobian[m];
		
	//delete []jacobian;
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Psi matrix related functions: Deleted because of the mistake, will be returned from the other program.

// Additional useful functions:

void Globals::set_w (double intw){
	
	w = intw;
}

double Globals::get_w (){
	
	return w;
}

void Globals::set_factor (int infactor){
	
	factor = infactor;
}

int Globals::get_factor (){
	
	return factor;
}
