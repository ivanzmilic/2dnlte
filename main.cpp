#include "all.h"

using namespace std;

int main() {
	
// This is the main function of the FBILI2D code. Here we will set-up the model and perform iterative solution of the RT problem.
// ------------------------------------------------------------------------------------------------------------------------------

// First we set up the geometry of the problem.

	int NX = 81;
	int NY = 81;
	Globals::setNX(NX);
	Globals::setNY(NY);
	
	int i, j, g, m; // Just in case

	Point ** Grid;
	Grid = new (nothrow) Point * [NX];
	for (i=0; i<NX; i++)
		Grid[i] = new (nothrow) Point [NY];
		
// -----------------------------------------------------------------------------------------------------------------------------
	// Set-up the geometry:
	
	if (geometry (Grid) == 0)
		cout << "Geometry set-up." << endl;
	
// -----------------------------------------------------------------------------------------------------------------------------

	// Set-up radiative transfer parameters:
	
	if (RT_parameters(Grid) == 0)
		cout << "Radiative transfer parameters set-up." << endl;
		
	int N = 21;
	Globals::allocate_x(N);
	Globals::compute_x();
	Globals::compute_profile();
	int n;
			
// ----------------------------------------------------------------------------------------------------------------------------

	// Here we setup the angles:
	
	Globals::compute_angles_Carlson ();
	//Globals::compute_angles_1();
	int NT = Globals::getNT();
	
// ----------------------------------------------------------------------------------------------------------------------------

	// Before the iteration begins, allocate I;
	
	for (i=0; i<NX; i++)
		for (j=0; j<NY; j++)
			Grid[i][j].allocate_I();

	int MAX_ITER = 500; // Maximum number of iterations we want to perform
	double sigma = 1E-4; // We stop the iterative process when we reach this relative change
	
// --------------------------------------------------------------------------------------------------------------------------------

	//////////////////  HERE THE ITERATIVE PART IS EXECUTED: //////////////////////////////////////////////////////////////////////
		
	
	Globals::set_w(1.0);
	Globals::set_factor(1.0);
	
	//Jacobi_explicit_full (Grid, MAX_ITER, sigma);	
	//GS_explicit_full (Grid, MAX_ITER, sigma);
	SSOR_explicit_full (Grid, 1000, sigma);
	//FBILI_explicit_full_2factors (Grid, MAX_ITER, sigma);
	
	//Globals::set_factor(1);
	
	
	ofstream output;
	output.open ("source_backup.txt");
	output.precision(10);
	for (i=0; i<NX; i++)
		for (j=0; j<NY; j++)
			output << i << " " << Grid[i][j].X << " " << j << " " << Grid[i][j].Y << " " << Grid[i][j].S << endl;
	output.close();

	output.open ("output_intensity.txt");
	output.precision(10);
	j = 0;
	cout << Globals::getNT() << " " << Globals::getNP() << " " << Globals::getN() << std::endl;
	for (i=0; i<NX; i++)
		for (g=0;g<Globals::getNT();++g)
				for (m=0;m<Globals::getNP_C(g);++m)
					for (int n=0;n<Globals::getN();++n)
						output << i << " " << Grid[i][j].X << " " << Grid[i][j].I[g][m][n] << endl;
	output.close();

		
// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	cout << "Execution finished." << endl;
	
	// Here we deallocate things:
	
	for (i=0; i<NX; i++)
		for (j=0; j<NY; j++)
			Grid[i][j].deallocate_I ();
	
	for (i = 0; i<NX; i++)
		delete []Grid[i];
		
	delete []Grid;
	
	cout << "Grid deallocated." << endl;
	
	return 0;
}
