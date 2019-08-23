#include "all.h"

using namespace std;

int main() {
	
// This is the main function of the FBILI2D code. Here we will set-up the model and perform iterative solution of the RT problem.
// ------------------------------------------------------------------------------------------------------------------------------

// First we set up the geometry of the problem.

	int NX = 123;
	int NY = 123;
	Globals::setNX(NX);
	Globals::setNY(NY);
	
	int i, j, g, m;

	Point ** Grid;
	Grid = new (nothrow) Point * [NX];
	for (i=0; i<NX; i++)
		Grid[i] = new (nothrow) Point [NY];
		
// -----------------------------------------------------------------------------------------------------------------------------
	// Set-up the geometry:
	
	if (geometry (Grid) == 0)
		cout << "Geometry set-up." << endl;
		
	//for (j=0; j<NY; j++)
		//cout << j << "\t" << Grid[0][j].Y << endl;
	
// -----------------------------------------------------------------------------------------------------------------------------

	// Set-up radiative transfer parameters:
	
	if (RT_parameters(Grid) == 0)
		cout << "Radiative transfer parameters set-up." << endl;
		
	int N = 9;
	Globals::allocate_x(N);
	Globals::compute_x();
	Globals::compute_profile();
	int n;
	//for (n=0; n<N; n++)
		//cout << n << "\t" << Globals::getx(n) << "\t" << Globals::getprofile(n) << "\t" << Globals::getwx(n) << endl;
	
		
// ----------------------------------------------------------------------------------------------------------------------------

	// Here we compute the angles:

	Globals::compute_angles_Carlson_10();	
	int NT = Globals::getNT();
	//cout << "NT = " << NT << endl;
	//for (g=0; g<NT; g++)
		//cout << Globals::getNP_C(g) << endl;
	
// ----------------------------------------------------------------------------------------------------------------------------

	// Before the iteration begins, allocate I;
	
	for (i=0; i<NX; i++)
		for (j=0; j<NY; j++)
			Grid[i][j].allocate_I();

	ofstream test;
	test.open ("Upwind_test.txt");
	n = 0;
	// Let us first carry out one iteration:
	
	double change;
	int i_max, j_max;
	int iter = 0;
	int MAX_ITER = 1000;
	
	int ray_counter = 0;
	
	int i_test = (NX-1)/2;
	int j_test = (NY-1)/2;
	
// --------------------------------------------------------------------------------------------------------------------------------

	//////////////////  HERE THE ITERATIVE PART IS EXECUTED: //////////////////////////////////////////////////////////////////////
	
	for (iter = 0; iter < MAX_ITER; iter ++){
		
		// Before anything else we compute the source function derivatives:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				Grid[i][j].compute_source_derivatives (Grid);
	
		// Direction ONE:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				// Allocate rays:
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=0; m<Globals::getNP_C(g)/4; m++){
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple (g, m, Grid);
						Grid[i][j].compute_pqr (g , m);
						Grid[i][j].formal_solution (g, m);
					}
				
				// Deallocate rays:
				Grid[i][j].deallocate_rays ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rxry ();
			}
			
		// Direction TWO:
		
		for (i= NX-1; i>-1; i--)
			for (j=0; j<NY; j++){
				
				// Allocate rays:
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g)/4; m<Globals::getNP_C(g)/2; m++){
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple (g, m, Grid);
						Grid[i][j].compute_pqr (g , m);
						Grid[i][j].formal_solution (g, m);
					}
				
				// Deallocate rays:
				Grid[i][j].deallocate_rays ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rxry ();
				
			}
			
		// Direction THREE:
		
		for (i= NX-1; i>-1; i--)
			for (j=NY-1; j>-1; j--){
				
				// Allocate rays:
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g)/2; m< 3 * Globals::getNP_C(g)/4; m++){
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple (g, m, Grid);
						Grid[i][j].compute_pqr (g , m);
						Grid[i][j].formal_solution (g, m);
					}
				
				// Deallocate rays:
				Grid[i][j].deallocate_rays ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rxry ();
				
			}
		
		
		// Direction FOUR:
			
		for (i= 0; i<NX; i++)
			for (j=NY-1; j>-1; j--){
	
				// Allocate rays:
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m= 3 * Globals::getNP_C(g)/4; m< Globals::getNP_C(g); m++){
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple (g, m, Grid);
						Grid[i][j].compute_pqr_return_variant1 (g , m);
						//Grid[i][j].compute_pqr_return_variant2 (g , m);
						//Grid[i][j].compute_pqr (g , m);
						//Grid[i][j].compute_pqr_linear (g , m );
						Grid[i][j].formal_solution_return_variant1 (g, m, Grid);
						//Grid[i][j].formal_solution_return_variant2 (g, m, Grid);
						//Grid[i][j].formal_solution (g, m);
					}
					
				Grid[i][j].compute_S_alo ();
				
				// Now updating:
				for (g=0; g<NT; g++)
					for (m= 3 * Globals::getNP_C(g)/4; m< Globals::getNP_C(g); m++){
						Grid[i][j].update_intensity (g, m);
					}
				
				// Deallocate rays:
				Grid[i][j].deallocate_rays ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rxry ();
				
			}
			
		// And then we can compute the new values of the source function:
		change = 0;
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				//Grid[i][j].compute_S_alo(); // Compute the source function
				if (aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS)) > aps (change)){
					change = Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS);
					i_max = i;
					j_max = j;
				}
			}
				
		cout << "Iteration number " << iter+1 << " finished. Maximum relative change is: " << change << " in the point : " << i_max << " , " << j_max << endl;
		
		if (aps(change) < 1E-5)
			break;
	}
	
	cout << "Iterative part finished. Total number of iterations is: " << iter << endl;

	test.close();
	
	ofstream output;
	output.open ("source.txt");
	for (i=0; i<NX; i++)
		for (j=0; j<NY; j++)
			output << i << " " << Grid[i][j].X << " " << j << " " << Grid[i][j].Y << " " << Grid[i][j].S << " " << Grid[i][j].BB << endl;
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
