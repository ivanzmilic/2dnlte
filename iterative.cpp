#include "all.h"

using namespace std;

double FBILI_2by2 (Point ** Grid, int N_iterations, double delta){

	int iter; // number of current iteration
	double change = 0; // 
	int imax, jmax;
	
	int m, g, i, j; // For going with the flow ;-)
	
	// Some global values we are going to need:
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int NT = Globals::getNT();
	
	// Start the iteration:
	
	for (iter = 0; iter < N_iterations; iter ++){
		
		// First the derivatives are computed:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				Grid[i][j].compute_source_derivatives_mod(Grid);
				
		// DIRECTION ONE:
		/*
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays();
				Grid[i][j].allocate_pqr();
				Grid[i][j].allocate_rxry();
				
				for (g=0; g<NT; g++)
					for (m=0; m<Globals::getNP_C(g)/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr(g, m);
						Grid[i][j].formal_solution(g, m);
						
					}
					
				Grid[i][j].deallocate_rxry ();
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
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) / 2; m<Globals::getNP_C(g) * 3/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr(g, m);
						Grid[i][j].formal_solution(g,m);
						
					}
					
				Grid[i][j].deallocate_rxry ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
			}
		*/	
		// And, for the sake of testing:
		
		// DIRECTION TWO:
			
		for (i=NX-1; i>-1; i--)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr(g, m);
						Grid[i][j].formal_solution(g,m);
						
					}
					
				Grid[i][j].deallocate_rxry ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
		
		// DIRECTION FOUR:	
			
		for (i=0; i<NX; i++)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_4th_direction(g, m);
						Grid[i][j].formal_solution(g,m);
						
					}
					
				// Now there is a need to compute modified derivatives in a different way and to execute return iterations. 
				
				Grid[i][j].modify_derivatives_asymmetric (Grid);
				
				Grid[i][j].compute_S_simple();
				
				// And to modify the local intensity:
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++)
						for (int n=0; n<Globals::getN(); n++)
							Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
				
				
				Grid[i][j].deallocate_rxry ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
				
		//for (i=0; i<NX; i++)
			//for (j=0; j<NY; j++)
				//Grid[i][j].compute_S_simple();
				
		change = 0.0;
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				if (aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS)) > change)
					change = aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS));
				
		cout << "Iteration number " << iter + 1 << " finished. Maximum relative change is :" << change << endl;
		
		if (change < delta)
			break;
				
	}
		
	return 0;
}

double FBILI_2by2_explicit_derivatives (Point ** Grid, int N_iterations, double delta){

	int iter; // number of current iteration
	double change = 0; // 
	int imax, jmax;
	
	int m, g, i, j; // For going with the flow ;-)
	
	// Some global values we are going to need:
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int NT = Globals::getNT();
	
	// Start the iteration:
	
	for (iter = 0; iter < N_iterations; iter ++){
		
		// First the derivatives are computed:
				
		// DIRECTION ONE:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays();
				Grid[i][j].allocate_pqr();
				Grid[i][j].allocate_rs();
				
				for (g=0; g<NT; g++)
					for (m=0; m<Globals::getNP_C(g)/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr_explicit_derivatives (g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
					
				if (iter > 0)
					Grid[i][j].compute_S_simple_explicit_derivatives (Grid);
				
				if (iter > 0){
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 0/4; m<Globals::getNP_C(g) * 1/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
							}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
				Grid[i][j].AA3 = 0.0;
			}
				
		// Now, contrary to the usual procedure, we shall sweep the grid in the direction 3, that is the direction opposite to the 
		// First sweep:
		
		// DIRECTION THREE:
		
		for (i=NX-1; i>-1; i--)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) / 2; m<Globals::getNP_C(g) * 3/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
				
				if (iter > 0)
					Grid[i][j].compute_S_simple_explicit_derivatives (Grid);
				
				if (iter > 0){
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 2/4; m<Globals::getNP_C(g) * 3/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
							}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
				Grid[i][j].AA2 = 0.0;
			}
			
		// And, for the sake of testing:
		
		// DIRECTION TWO:
			
		for (i=NX-1; i>-1; i--)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g,m, Grid, iter);
						
					}
					
				if (iter > 0)
					Grid[i][j].compute_S_simple_explicit_derivatives (Grid);
				
				if (iter > 0){
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
							}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
				Grid[i][j].AA4 = 0.0;
				
			}
		
		// DIRECTION FOUR:	
			
		for (i=0; i<NX; i++)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
					
				// Now there is a need to compute modified derivatives in a different way and to execute return iterations. 
				
				
				Grid[i][j].compute_S_simple_explicit_derivatives(Grid);
				
				// And to modify the local intensity:
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++)
						for (int n=0; n<Globals::getN(); n++)
							Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
				
				
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
			
		//for (i=0; i<NX; i++)
			//for (j=0; j<NY; j++)
				//Grid[i][j].compute_S_simple_explicit_derivatives(Grid);
				
		change = 0.0;
		imax = 0;
		jmax = 0;
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				if (aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS)) > change){
					change = aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS));
					imax = i;
					jmax = j;
				}
				
		cout << "Iteration number " << iter + 1 << " finished. Maximum relative change is :" << change << " in point: " 
			<< imax << " , " << jmax << endl;
		
		if (change < delta)
			break;
			
		// And here we reset the values of the coefficients:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				//Grid[i][j].BB = 0.0;
				
				Grid[i][j].AA1 = 0.0;
				//Grid[i][j].AA2 = 0.0;
				//Grid[i][j].AA3 = 0.0;
				//Grid[i][j].AA4 = 0.0;
				
				//Grid[i][j].C1 = 0.0;
				//Grid[i][j].C2 = 0.0;
				//Grid[i][j].C3 = 0.0;
				//Grid[i][j].C4 = 0.0;
			}
				
	}
		
	return 0;
}

double SSOR_explicit_derivatives (Point ** Grid, int N_iterations, double delta){

	int iter; // number of current iteration
	double change = 0; // 
	int imax, jmax;
	
	int m, g, i, j; // For going with the flow ;-)
	
	// Some global values we are going to need:
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int NT = Globals::getNT();
	
	// Start the iteration:
	
	for (iter = 0; iter < N_iterations; iter ++){
		
		// First the derivatives are computed:
				
		// DIRECTION ONE:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays();
				Grid[i][j].allocate_pqr();
				Grid[i][j].allocate_rs();
				
				for (g=0; g<NT; g++)
					for (m=0; m<Globals::getNP_C(g)/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr_explicit_derivatives (g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
					
				if (iter > 0)
					Grid[i][j].compute_S_simple_explicit_derivatives (Grid);
				
				if (iter > 0){
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 0/4; m<Globals::getNP_C(g) * 1/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
							}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
				Grid[i][j].AA3 = 0.0;
			}
				
		// Now, contrary to the usual procedure, we shall sweep the grid in the direction 3, that is the direction opposite to the 
		// First sweep:
		
		// DIRECTION THREE:
		
		for (i=NX-1; i>-1; i--)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) / 2; m<Globals::getNP_C(g) * 3/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
				
				if (iter > 0)
					Grid[i][j].compute_S_simple_explicit_derivatives (Grid);
				
				if (iter > 0){
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 2/4; m<Globals::getNP_C(g) * 3/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
							}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
				Grid[i][j].AA2 = 0.0;
			}
			
		// And, for the sake of testing:
		
		// DIRECTION TWO:
			
		for (i=NX-1; i>-1; i--)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g,m, Grid, iter);
						
					}
					
				if (iter > 0)
					Grid[i][j].compute_S_simple_explicit_derivatives (Grid);
				
				if (iter > 0){
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
							}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
				Grid[i][j].AA4 = 0.0;
				
			}
		
		// DIRECTION FOUR:	
			
		for (i=0; i<NX; i++)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
					
				// Now compute new values for the Source function:
				
				Grid[i][j].compute_S_simple_explicit_derivatives(Grid);
				
				// And to modify the local intensity:
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++)
						for (int n=0; n<Globals::getN(); n++)
							Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
				
				
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
			
		//for (i=0; i<NX; i++)
			//for (j=0; j<NY; j++)
				//Grid[i][j].compute_S_simple_explicit_derivatives(Grid);
				
		change = 0.0;
		imax = 0;
		jmax = 0;
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				if (aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS)) > change){
					change = aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS));
					imax = i;
					jmax = j;
				}
				
		cout << "Iteration number " << iter + 1 << " finished. Maximum relative change is :" << change << " in point: " 
			<< imax << " , " << jmax << endl;
		
		if (change < delta)
			break;
			
		// And here we reset the values of the coefficients:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				//Grid[i][j].BB = 0.0;
				
				Grid[i][j].AA1 = 0.0;
				//Grid[i][j].AA2 = 0.0;
				//Grid[i][j].AA3 = 0.0;
				//Grid[i][j].AA4 = 0.0;
				
				//Grid[i][j].C1 = 0.0;
				//Grid[i][j].C2 = 0.0;
				//Grid[i][j].C3 = 0.0;
				//Grid[i][j].C4 = 0.0;
			}
				
	}
		
	return 0;
}

double Jacobi_explicit_derivatives (Point ** Grid, int N_iterations, double delta){

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
	convergence.open("convergence_jacobi.txt");
	
	for (iter = 0; iter < N_iterations; iter ++){
		
		// First the derivatives are computed:
				
		// DIRECTION ONE:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays();
				Grid[i][j].allocate_pqr();
				Grid[i][j].allocate_rs();
				
				for (g=0; g<NT; g++)
					for (m=0; m<Globals::getNP_C(g)/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr_explicit_derivatives (g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
						
				Grid[i][j].deallocate_rs ();
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
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) / 2; m<Globals::getNP_C(g) * 3/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
			
		// And, for the sake of testing:
		
		// DIRECTION TWO:
			
		for (i=NX-1; i>-1; i--)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g,m, Grid, iter);
						
					}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
		
		// DIRECTION FOUR:	
			
		for (i=0; i<NX; i++)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
				
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
			
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				Grid[i][j].add_a_and_c(Grid);
			
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				Grid[i][j].compute_S_simple();
				
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
		convergence << iter << " " << change << " " << Grid[(NX-1)/2][0].S << endl;
				
		cout << "Iteration number " << iter + 1 << " finished. Maximum relative change is :" << change << " in point: " 
			<< imax << " , " << jmax << endl;
		
		if (aps(change) < delta)
			break;
			
		// And here we reset the values of the coefficients:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].AA1 = 0.0;
				Grid[i][j].AA2 = 0.0;
				Grid[i][j].AA3 = 0.0;
				Grid[i][j].AA4 = 0.0;
				
			}
				
	}
	
	convergence.close();
		
	return 0;
}

double GS_explicit_derivatives (Point ** Grid, int N_iterations, double delta){

	int iter; // number of current iteration
	double change = 0; // 
	int imax, jmax;
	
	int m, g, i, j; // For going with the flow ;-)
	
	// Some global values we are going to need:
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int NT = Globals::getNT();
	
	// Start the iteration:
	
	for (iter = 0; iter < N_iterations; iter ++){
		
		// First the derivatives are computed:
				
		// DIRECTION ONE:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays();
				Grid[i][j].allocate_pqr();
				Grid[i][j].allocate_rs();
				
				for (g=0; g<NT; g++)
					for (m=0; m<Globals::getNP_C(g)/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr_explicit_derivatives (g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
						
				Grid[i][j].deallocate_rs ();
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
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) / 2; m<Globals::getNP_C(g) * 3/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
			
		// And, for the sake of testing:
		
		// DIRECTION TWO:
			
		for (i=NX-1; i>-1; i--)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g,m, Grid, iter);
						
					}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
		
		// DIRECTION FOUR:	
			
		for (i=0; i<NX; i++)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives(g, m, Grid, iter);
						
					}
					
				// Now compute new values for the Source function:
				
				Grid[i][j].compute_S_simple_explicit_derivatives(Grid);
				
				// And to modify the local intensity:
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++)
						for (int n=0; n<Globals::getN(); n++)
							Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
				
				
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
			
		//for (i=0; i<NX; i++)
			//for (j=0; j<NY; j++)
				//Grid[i][j].compute_S_simple_explicit_derivatives(Grid);
				
		change = 0.0;
		imax = 0;
		jmax = 0;
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				if (aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS)) > change){
					change = aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS));
					imax = i;
					jmax = j;
				}
				
		cout << "Iteration number " << iter + 1 << " finished. Maximum relative change is :" << change << " in point: " 
			<< imax << " , " << jmax << endl;
		
		if (change < delta)
			break;
			
		// And here we reset the values of the coefficients:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				//Grid[i][j].BB = 0.0;
				
				Grid[i][j].AA1 = 0.0;
				Grid[i][j].AA2 = 0.0;
				Grid[i][j].AA3 = 0.0;
				Grid[i][j].AA4 = 0.0;
				
				//Grid[i][j].C1 = 0.0;
				//Grid[i][j].C2 = 0.0;
				//Grid[i][j].C3 = 0.0;
				//Grid[i][j].C4 = 0.0;
			}
				
	}
		
	return 0;
}

double FBILI_explicit_derivatives (Point ** Grid, int N_iterations, double delta){

	int iter; // number of current iteration
	double change = 0; // 
	int imax, jmax;
	
	int m, g, i, j; // For going with the flow ;-)
	
	// Some global values we are going to need:
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int NT = Globals::getNT();
	
	// Start the iteration:
	
	for (iter = 0; iter < N_iterations; iter ++){
		
		// First the derivatives are computed:
				
		// DIRECTION ONE:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays();
				Grid[i][j].allocate_pqr();
				Grid[i][j].allocate_rs();
				
				for (g=0; g<NT; g++)
					for (m=0; m<Globals::getNP_C(g)/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr_explicit_derivatives (g, m);
						Grid[i][j].formal_solution_explicit_derivatives_factors(g, m, Grid, iter);
						
					}
						
				Grid[i][j].deallocate_rs ();
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
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) / 2; m<Globals::getNP_C(g) * 3/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives_factors(g, m, Grid, iter);
						
					}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
			
		// And, for the sake of testing:
		
		// DIRECTION TWO:
			
		for (i=NX-1; i>-1; i--)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives_factors(g,m, Grid, iter);
						
					}
					
				Grid[i][j].deallocate_rs ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
		
		// DIRECTION FOUR:	
			
		for (i=0; i<NX; i++)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rs ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_explicit_derivatives(g, m);
						Grid[i][j].formal_solution_explicit_derivatives_factors(g, m, Grid, iter);
						
					}
					
				// Now compute new values for the Source function:
				
				Grid[i][j].compute_S_simple_explicit_derivatives_factors(Grid);
				
				// And to modify the local intensity:
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++)
						for (int n=0; n<Globals::getN(); n++)
							Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
				
				
				Grid[i][j].deallocate_rs ();
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
		
		if (aps(change) < delta)
			break;
			
		// And here we reset the values of the coefficients:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				//Grid[i][j].BB = 0.0;
				
				Grid[i][j].AA1 = 0.0;
				Grid[i][j].AA2 = 0.0;
				Grid[i][j].AA3 = 0.0;
				Grid[i][j].AA4 = 0.0;
				
				Grid[i][j].BB1 = 0.0;
				Grid[i][j].BB2 = 0.0;
				Grid[i][j].BB3 = 0.0;
				Grid[i][j].BB4 = 0.0;
						
				//Grid[i][j].C1 = 0.0;
				//Grid[i][j].C2 = 0.0;
				//Grid[i][j].C3 = 0.0;
				//Grid[i][j].C4 = 0.0;
			}
				
	}
		
	return 0;
}

// Fully explicit versions:

double Jacobi_explicit_full (Point ** Grid, int N_iterations, double delta){

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
	convergence.open("convergence_J.txt");
	
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
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);
						
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);
						
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g,m, Grid, iter);
						
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);
						
					}
					
				// Now compute new values for the Source function:
				
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}

		for (i=0; i<NX; i++)
			for (j=0; j<NY;j++)
				Grid[i][j].modify_a(Grid);
			
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				Grid[i][j].compute_S_simple();
				
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

		ofstream debug;
		debug.open("debug.txt");

		/*for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				debug << i << " " << j << " " << Grid[i][j].AA1 << " " << Grid[i][j].AA2 << " " << Grid[i][j].AA3 << " " << Grid[i][j].AA4 << endl;
				debug << i << " " << j << " " << Grid[i][j].BB1 << " " << Grid[i][j].BB2 << " " << Grid[i][j].BB3 << " " << Grid[i][j].BB4 << endl; 
				debug << i << " " << j << " " << Grid[i][j].C1 << " " << Grid[i][j].C2 << " " << Grid[i][j].C3 << " " << Grid[i][j].C4 << endl; 
				debug << i << " " << j << " " << Grid[i][j].C5 << " " << Grid[i][j].C6 << " " << Grid[i][j].C7 << " " << Grid[i][j].C8 << endl; 
			}*/

		debug.close();
		
		if (aps(change) < delta)
			break;
			
		// And here we reset the values of the coefficients:

		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].BB1 = 0.0;
				Grid[i][j].BB2 = 0.0;
				Grid[i][j].BB3 = 0.0;
				Grid[i][j].BB4 = 0.0;
				
				Grid[i][j].AA1 = 0.0;
				Grid[i][j].AA2 = 0.0;
				Grid[i][j].AA3 = 0.0;
				Grid[i][j].AA4 = 0.0;
				
				//Grid[i][j].C1 = 0.0;
				//Grid[i][j].C2 = 0.0;
				//Grid[i][j].C3 = 0.0;
				//Grid[i][j].C4 = 0.0;
			}
				
	}
	
	convergence.close();
		
	return 0;
}

double GS_explicit_full (Point ** Grid, int N_iterations, double delta){

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
	convergence.open("convergence_GS.txt");
	
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
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);
						
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);
						
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g,m, Grid, iter);
						
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);
						
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
			
		//for (i=0; i<NX; i++)
			//for (j=0; j<NY; j++)
				//Grid[i][j].compute_S_simple_explicit_full(Grid);
				
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
				Grid[i][j].BB2 = 0.0;
				Grid[i][j].BB3 = 0.0;
				Grid[i][j].BB4 = 0.0;
				
				Grid[i][j].AA1 = 0.0;
				Grid[i][j].AA2 = 0.0;
				Grid[i][j].AA3 = 0.0;
				Grid[i][j].AA4 = 0.0;
				
				//Grid[i][j].C1 = 0.0;
				//Grid[i][j].C2 = 0.0;
				//Grid[i][j].C3 = 0.0;
				//Grid[i][j].C4 = 0.0;
			}
				
	}
	
	convergence.close();
		
	return 0;
}

double SSOR_explicit_full (Point ** Grid, int N_iterations, double delta){

	int iter; // number of current iteration
	double change = 0; // 
	int imax, jmax;
	
	int m, g, i, j; // For going with the flow ;-)
	
	// Some global values we are going to need:
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int NT = Globals::getNT();

	int limit = 0;
	
	// Start the iteration:
	
	ofstream convergence;
	convergence.precision(10);
	convergence.open("convergence_SSOR.txt");
	
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
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);		
					}
					
				if (iter > limit){
					// Now compute new values for the Source function:
		
					Grid[i][j].compute_S_simple_explicit_full(Grid);
					
					// And to modify the local intensity:
					
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 0/4; m<Globals::getNP_C(g) * 1/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
					
					
				}
						
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
				if (iter >0){
					Grid[i][j].AA3 = 0.0;
					Grid[i][j].BB3 = 0.0;
				}
				
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);
						
					}
					
				if (iter > limit){
					// Now compute new values for the Source function:
		
					Grid[i][j].compute_S_simple_explicit_full(Grid);
					
					// And to modify the local intensity:
					
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 2/4; m<Globals::getNP_C(g) * 3/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
						
				}
					
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
				if (iter >0){
					Grid[i][j].AA2 = 0.0;
					Grid[i][j].BB2 = 0.0;
				}
				
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g,m, Grid, iter);
						
					}
					
				if (iter > limit){
					// Now compute new values for the Source function:
		
					Grid[i][j].compute_S_simple_explicit_full(Grid);
					
					// And to modify the local intensity:
					
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;	
				}
					
				Grid[i][j].deallocate_r_full ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
				// And reset appropriate coefficients:
				if (iter >0){
					Grid[i][j].AA4 = 0.0;
					Grid[i][j].BB4 = 0.0;
				}
				
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);
						
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
				
				// And reset appropriate coefficients:
				
				Grid[i][j].AA1 = 0.0;
				Grid[i][j].BB1 = 0.0;
				
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
			
	}
		
	convergence.close();
	
	return 0;
}

double FBILI_explicit_full_2factors (Point ** Grid, int N_iterations, double delta){

	int iter; // number of current iteration
	double change = 0; // 
	int imax, jmax;
	
	int m, g, i, j; // For going with the flow ;-)
	
	// Some global values we are going to need:
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int NT = Globals::getNT();

	int limit = 0;
	
	// Start the iteration:
	
	ofstream convergence;
	convergence.precision(10);
	convergence.open("convergence_FBILI.txt");
	ofstream twodplots;
	twodplots.precision(7);


	Globals::set_factor(1);
	
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
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);
						
					}
				
				if (iter > limit){
					
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);
						
					}
					
				if (iter > limit){
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g,m, Grid, iter);
						
					}
				
				if (iter > limit){
					
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full(g, m, Grid, iter);
						
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
			
		//for (i=0; i<NX; i++)
			//for (j=0; j<NY; j++)
				//Grid[i][j].compute_S_simple_explicit_full(Grid);
				
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



		if (iter == 9){
			twodplots.open("fbili_10.txt");
			for (i=0; i<NX; i++)
				for (j=0; j<NY; j++)
					twodplots << i << " " << Globals::getX(i) << " " << j << " " << Globals::getY(j) << " " << Grid[i][j].S << endl;
			twodplots.close();
		}

		if (iter == 19){
			twodplots.open("fbili_20.txt");
			for (i=0; i<NX; i++)
				for (j=0; j<NY; j++)
					twodplots << i << " " << Globals::getX(i) << " " << j << " " << Globals::getY(j) << " " << Grid[i][j].S << endl;
			twodplots.close();
		}
		
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
	
	for (i=0; i<NX; i++)
		for (j=0; j<NY; j++){
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
		}
	
	return 0;
}

double FBILI_explicit_full_2by2 (Point ** Grid, int N_iterations, double delta){

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
	convergence.open("convergence_FBILI_2by2.txt");
	
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
						Grid[i][j].formal_solution_explicit_full_factors(g, m, Grid, iter);
						
					}
				
				if (iter > 0){
					
					/*Grid[i][j].compute_S_simple_explicit_full(Grid);
					
					// And to modify the local intensity:
					
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 0/4; m<Globals::getNP_C(g) * 1/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;*/
					
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full_factors(g, m, Grid, iter);
						
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full_factors(g,m, Grid, iter);
						
					}
				
				if (iter > 0){
					
					/*Grid[i][j].compute_S_simple_explicit_full(Grid);
					
					// And to modify the local intensity:
					
					for (g=0; g<NT; g++)
						for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++)
							for (int n=0; n<Globals::getN(); n++)
								Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;*/
						
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
						Grid[i][j].compute_pqr_explicit_full(g, m);
						Grid[i][j].formal_solution_explicit_full_factors(g, m, Grid, iter);
						
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
			
		//for (i=0; i<NX; i++)
			//for (j=0; j<NY; j++)
				//Grid[i][j].compute_S_simple_explicit_full(Grid);
				
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
	
	for (i=0; i<NX; i++)
		for (j=0; j<NY; j++){
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
		}
	
	return 0;
}

// ---------------------------------------------------------------------------------------------

// "Old style"

double FBILI_upwind_derivative (Point ** Grid, int N_iterations, double delta){
	
	int iter; // number of current iteration
	double change = 0; // 
	int imax, jmax;
	
	int m, g, i, j; // For going with the flow ;-)
	
	// Some global values we are going to need:
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int NT = Globals::getNT();
	
	// Start the iteration:
	
	for (iter = 0; iter < 1; iter ++){
		
		// First the derivatives are computed:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				Grid[i][j].compute_source_derivatives_mod(Grid);
				
		// DIRECTION ONE:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays();
				Grid[i][j].allocate_pqr();
				Grid[i][j].allocate_rxry();
				
				for (g=0; g<NT; g++)
					for (m=0; m<Globals::getNP_C(g)/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr(g, m);
						Grid[i][j].formal_solution(g, m);
						
					}
					
				Grid[i][j].deallocate_rxry ();
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
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) / 2; m<Globals::getNP_C(g) * 3/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr(g, m);
						Grid[i][j].formal_solution(g,m);
						
					}
					
				Grid[i][j].deallocate_rxry ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
			}

		// And, for the sake of testing:
		
		// DIRECTION TWO:
			
		for (i=NX-1; i>-1; i--)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr(g, m);
						Grid[i][j].formal_solution(g,m);
						
					}
					
				Grid[i][j].deallocate_rxry ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
		
		// DIRECTION FOUR:	
			
		for (i=0; i<NX; i++)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m, Grid);
						Grid[i][j].compute_pqr_4th_direction(g, m);
						Grid[i][j].formal_solution(g,m);
						
					}
					
				// Now there is a need to compute modified derivatives in a different way and to execute return iterations. 
				
				Grid[i][j].modify_derivatives_asymmetric (Grid);
				
				Grid[i][j].compute_S_simple();
				
				// And to modify the local intensity:
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++)
						for (int n=0; n<Globals::getN(); n++)
							Grid[i][j].I[g][m][n] += Grid[i][j].p[g][m][n] * Grid[i][j].dS;
				
				
				Grid[i][j].deallocate_rxry ();
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				
			}
				
				
		change = 0.0;
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				if (aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS)) > change)
					change = aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS));
				
		cout << "Iteration number " << iter + 1 << " finished. Maximum relative change is :" << change << endl;
		
		if (change < delta)
			break;
				
	}
		
	return 0;	
}


// Additional functions to read S from file.

int read_S_from_file (Point ** Grid){
	
	ifstream input("source_backup.txt",std::ios::in |std::ios::binary);
	
	int i, j;
	
	int NX = Globals::getNX();
	int NY = Globals::getNY();
	
	double temp;
	
	for (i=0; i<NX; i++)
		for (j=0; j<NY; j++){
			input >> temp;
			input >> temp;
			input >> temp;
			input >> temp;
			input >> temp;
			
			Grid[i][j].S = temp;
		}
			
	return 0;
}

// Re-attempt @ elmination scheme

double Jacobi_classic (Point ** Grid, int N_iterations, double delta){
	
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
	convergence.open("convergence_test.txt");
	
	for (iter = 0; iter < N_iterations; iter ++){
		
		// First the derivatives are computed:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++)
				Grid[i][j].compute_source_derivatives(Grid);
				
		// DIRECTION ONE:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays();
				Grid[i][j].allocate_pqr();
				Grid[i][j].allocate_rxry();
				
				for (g=0; g<NT; g++)
					for (m=0; m<Globals::getNP_C(g)/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr_classic(g, m);
						Grid[i][j].formal_solution_classic(g, m, Grid, iter);
						
					}
						
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				Grid[i][j].deallocate_rxry ();
				
			}
				
		// Now, contrary to the usual procedure, we shall sweep the grid in the direction 3, that is the direction opposite to the 
		// First sweep:
		
		// DIRECTION THREE:
		
		for (i=NX-1; i>-1; i--)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) / 2; m<Globals::getNP_C(g) * 3/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr_classic(g, m);
						Grid[i][j].formal_solution_classic(g, m, Grid, iter);
						
					}
					
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				Grid[i][j].deallocate_rxry ();
				
			}
			
		// DIRECTION TWO:
			
		for (i=NX-1; i>-1; i--)
			for (j=0; j<NY; j++){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 1/4; m<Globals::getNP_C(g) * 2/4; m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr_classic(g, m);
						Grid[i][j].formal_solution_classic(g, m, Grid, iter);
						
					}
					
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				Grid[i][j].deallocate_rxry ();
				
			}
		
		// DIRECTION FOUR:	
			
		for (i=0; i<NX; i++)
			for (j=NY-1; j>-1; j--){
				
				Grid[i][j].allocate_rays ();
				Grid[i][j].allocate_pqr ();
				Grid[i][j].allocate_rxry ();
				
				for (g=0; g<NT; g++)
					for (m=Globals::getNP_C(g) * 3/4; m<Globals::getNP_C(g); m++){
						
						Grid[i][j].compute_upwind_points(g, m, Grid);
						Grid[i][j].interpolate_simple(g, m , Grid);
						Grid[i][j].compute_pqr_classic(g, m);
						Grid[i][j].formal_solution_classic(g, m, Grid, iter);
						
					}
				Grid[i][j].compute_S_classic(Grid);
				
				Grid[i][j].deallocate_pqr ();
				Grid[i][j].deallocate_rays ();
				Grid[i][j].deallocate_rxry ();
				
			}
			
			
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				//Grid[i][j].compute_S_classic(Grid);
			}
				
		change = 0.0;
		imax = 0;
		jmax = 0;
		
		for (i=59; i<60; i++)
			for (j=0; j<NY; j++)
				if (aps(Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS)) > aps(change)){
					change = (Grid[i][j].dS / (Grid[i][j].S - Grid[i][j].dS));
					imax = i;
					jmax = j;
				}
		convergence << iter << " " << change << " " << Grid[(NX-1)/2][0].S << endl;
				
		cout << "Iteration number " << iter + 1 << " finished. Maximum relative change is :" << change << " in point: " 
			<< imax << " , " << jmax << endl;
		
		if (aps(change) < delta)
			break;
			
		// And here we reset the values of the coefficients:
		
		for (i=0; i<NX; i++)
			for (j=0; j<NY; j++){
				
				Grid[i][j].AA = 0.0;
				Grid[i][j].BB = 0.0;
				Grid[i][j].CX = 0.0;
				Grid[i][j].CY = 0.0;
				
			}
				
	}
	
	convergence.close();
		
	return 0;
}


