#include "all.h"

using namespace std;

double FBILI_explicit_full_variant1 (Point ** Grid, int N_iterations, double delta){

	int iter; // number of current iteration
	double change = 0; // 
	int imax, jmax;
	
	int m, g, i, j; // For going with the flow ;-)
	
	// Some global values we are going to need:
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int NT = Globals::getNT();
	
	// Start the iteration:
	
	ofstream convergence;
	convergence.precision(10);
	convergence.open("convergence_FBILI_v1.txt");
	
	for (iter = 0; iter < N_iterations; iter ++){
		
		// First the derivatives are computed:
				
		// DIRECTION ONE:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays();
				Grid[i][j].allocate_pqr();
				Grid[i][j].allocate_r_full();
				
				for (g=0; g<NT; g++)
					for (m=0; m<Globals::getNP_C(g)/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr_explicit_full_v1 (g, m);
						Grid[i][j].formal_solution_explicit_full_v1(g, m, Grid, iter);
						
					}
				
				if (iter > 0){
					
					Grid[i][j].compute_S_simple_explicit_full(Grid);
					
					// And to modify the local intensity:
					
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 0/4; m<Globals::getNP_C(g) * 1/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
					
					Grid[i][j].AA3 = 0.0;
					Grid[i][j].BB3 = 0.0;
				}
						
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
				
		// Now, contrary to the usual procedure, we shall sweep the grid in the direction 3, that is the direction opposite to the 
		// First sweep:
		
		// DIRECTION THREE:
		
		for (i=NX-1; i>-1; i--)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_r_full ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) / 2; m<Globals::getNP_C(g) * 3/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_full_v1 (g, m);
						Grid[i][j].formal_solution_explicit_full_v1(g, m, Grid, iter);
						
					}
					
				if (iter > 0){
				// Now compute new values for the Source function:
	
					Grid[i][j].compute_S_simple_explicit_full(Grid);
					
					// And to modify the local intensity:
					
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 2/4; m<Globals::getNP_C(g) * 3/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
								
					Grid[i][j].AA2 = 0.0;
					Grid[i][j].BB2 = 0.0;
				}
					
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
			
		// And, for the sake of testing:
		
		// DIRECTION TWO:
			
		for (i=NX-1; i>-1; i--)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_r_full ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_full_v1 (g, m);
						Grid[i][j].formal_solution_explicit_full_v1(g, m, Grid, iter);
						
					}
				
				if (iter > 0){
					
					Grid[i][j].compute_S_simple_explicit_full(Grid);
					
					// And to modify the local intensity:
					
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
						
					Grid[i][j].AA4 = 0.0;
					Grid[i][j].BB4 = 0.0;
				}
					
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
		
		// DIRECTION FOUR:	
			
		for (i=0; i<NX; i++)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_r_full ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_full_v1 (g, m);
						Grid[i][j].formal_solution_explicit_full_v1(g, m, Grid, iter);
						
					}
					
				// Now compute new values for the Source function:
				
				Grid[i][j].compute_S_simple_explicit_full(Grid);
				
				// And to modify the local intensity:
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++)
						for (int n=0; n<Globals::getN(); n++)
							Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
				
				
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
			
	
		change = 0.0;
		imax = 0;
		jmax = 0;
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				if (aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS)) > aps(change)){
					change = (Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS));
					imax = i;
					jmax = j;
				}
				
		cout << "Iteration number " << iter + 1 << " finished. Maximum relative change is :" << change << " in point: " 
			<< imax << " , " << jmax << endl;
			
		convergence << iter << " " << change << " " << Grid[(NX-1)/2][0].S << endl;
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].BB1 = 0.0;
				
				Grid[i][j].AA1 = 0.0;
				
			}
		
		if (aps(change) < delta)
			break;
			
		// And here we reset the values of the coefficients:	
	}
		
	convergence.close();
	
	// Reset all the coefficients to zero in case you want to continue with another iterative scheme:
	
	return 0;
}

double FBILI_explicit_full_variant2 (Point ** Grid, int N_iterations, double delta){

	int iter; // number of current iteration
	double change = 0; // 
	int imax, jmax;
	
	int m, g, i, j; // For going with the flow ;-)
	
	// Some global values we are going to need:
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int NT = Globals::getNT();
	
	// Start the iteration:
	
	ofstream convergence;
	convergence.precision(10);
	convergence.open("convergence_FBILI_v1.txt");
	
	for (iter = 0; iter < N_iterations; iter ++){
		
		// First the derivatives are computed:
				
		// DIRECTION ONE:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays();
				Grid[i][j].allocate_pqr();
				Grid[i][j].allocate_r_full();
				
				for (g=0; g<NT; g++)
					for (m=0; m<Globals::getNP_C(g)/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr_explicit_full (g, m);
						Grid[i][j].formal_solution_explicit_full_v2(g, m, Grid, iter);
						
					}
				
				if (iter > -1){
					
					Grid[i][j].compute_S_simple_explicit_full(Grid);
					
					// And to modify the local intensity:
					
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 0/4; m<Globals::getNP_C(g) * 1/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
					
					Grid[i][j].AA3 = 0.0;
					Grid[i][j].BB3 = 0.0;
				}
						
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
		
		// Now, contrary to the usual procedure, we shall sweep the grid in the direction 3, that is the direction opposite to the 
		// First sweep:
		
		// DIRECTION THREE:
		
		for (i=NX-1; i>-1; i--)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_r_full ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) / 2; m<Globals::getNP_C(g) * 3/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_full (g, m);
						Grid[i][j].formal_solution_explicit_full_v2(g, m, Grid, iter);
						
					}
					
				if (iter > -1){
				// Now compute new values for the Source function:
	
					Grid[i][j].compute_S_simple_explicit_full(Grid);
					
					// And to modify the local intensity:
					
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 2/4; m<Globals::getNP_C(g) * 3/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
								
					Grid[i][j].AA2 = 0.0;
					Grid[i][j].BB2 = 0.0;
				}
					
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
			
		// And, for the sake of testing:
		
		// DIRECTION TWO:
			
		for (i=NX-1; i>-1; i--)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_r_full ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_full (g, m);
						Grid[i][j].formal_solution_explicit_full_v2(g, m, Grid, iter);
						
					}
				
				if (iter > -1){
					
					Grid[i][j].compute_S_simple_explicit_full(Grid);
					
					// And to modify the local intensity:
					
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
						
					Grid[i][j].AA4 = 0.0;
					Grid[i][j].BB4 = 0.0;
				}
					
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
		
		// DIRECTION FOUR:	
			
		for (i=0; i<NX; i++)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_r_full ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_full (g, m);
						Grid[i][j].formal_solution_explicit_full_v2(g, m, Grid, iter);
						
					}
					
				// Now compute new values for the Source function:
				
				Grid[i][j].compute_S_simple_explicit_full(Grid);
				
				// And to modify the local intensity:
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++)
						for (int n=0; n<Globals::getN(); n++)
							Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
				
				
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
							
		change = 0.0;
		imax = 0;
		jmax = 0;
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				if (aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS)) > aps(change)){
					change = (Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS));
					imax = i;
					jmax = j;
				}
				
		cout << "Iteration number " << iter + 1 << " finished. Maximum relative change is :" << change << " in point: " 
			<< imax << " , " << jmax << endl;
			
		convergence << iter << " " << change << " " << Grid[(NX-1)/2][0].S << endl;
		
		if (aps(change) < delta)
			break;
			
		// And here we reset the values of the coefficients:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].BB1 = 0.0;
				
				Grid[i][j].AA1 = 0.0;
				
			}
				
	}
		
	convergence.close();
	
	// Reset all the coefficients to zero in case you want to continue with another iterative scheme:
	
	for (i=0; i<NX; i++)
		for (j=0; j<NY; j++){
			/*
			Grid[i][j].AA1 = 0.0;
			Grid[i][j].AA2 = 0.0;
			Grid[i][j].AA3 = 0.0;
			Grid[i][j].AA4 = 0.0;
			
			Grid[i][j].BB1 = 0.0;
			Grid[i][j].BB2 = 0.0;
			Grid[i][j].BB3 = 0.0;
			Grid[i][j].BB4 = 0.0;
			
			Grid[i][j].C1 = 0.0;
			Grid[i][j].C2 = 0.0;
			Grid[i][j].C3 = 0.0;
			Grid[i][j].C4 = 0.0;
			Grid[i][j].C5 = 0.0;
			Grid[i][j].C6 = 0.0;
			Grid[i][j].C7 = 0.0;
			Grid[i][j].C8 = 0.0;
			*/
		}
	
	return 0;
}

