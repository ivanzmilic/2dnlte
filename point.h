#ifndef POINT_H
#define POINT_H

class Ray{
	
	// This is a class which holds all the needed interpolated quantities.
	
	private:
	
	// Nothing.
	
	public:
	
	// Geometrical coordinates of the upwind point:
	double Xu;
	double Yu;
	
	int type_u;
	double * dt_u; // geometrical path distance between local and the upwind point.
	double Chi; // opacity at the point
	
	// Physical quantities in the upwind point:
	
	double Su;
	double Spu;
	double * Iu;
	// Strangely I CANNOT remember if I need anything else but this :)
	
	// Procedures:
	int clear_ray ();
};

class Point{
private:
public:
	// Geometrical:
	double X; 
	double Y;
	int i, j;

	// Physical:

	double Chi; // Later computed from n
	double B; // Later replaced by/ computed from T
	double ep; // Later computed from n

	// Radiative transfer related:

	double S; // Source function, later emissivity
	double dS; // Correction for source function.
	double J; // Mean intensity, later replaced by irreducibile tensors.
	double alo; 
	double Spx; // Derivative along X
	double Spy; // Derivative along Y
	
	double Spx_mod; // Modified derivative along X (without the contribution of the local point)
	double Spy_mod; // Modified derivative along Y (without the contribution of the local point)

	double *** I;
	double *** p, *** q, *** r;
	double *** exxp;
	double *** rx, *** ry;
	double *** r1, *** r2, *** r3, *** r4;
	
	double *** r5, *** r6, *** r7, *** r8;
	
	// Additional variables for the polarization:
		
	double AA, BB;
	double CX, CY;
	double C1, C2, C3, C4;
	
	double C5, C6, C7, C8; // These are the additional coefficients which arise because we want to explicitly treat upwind Source function.
	
	double AA1, AA2, AA3, AA4; // AA broken in fragments, one for each sweep
	
	double BB1, BB2, BB3, BB4;
	
	// And, of course, rays! 
	
	Ray ** rays;

// ---------------------------------------------------------------------------------------------------------------------------------------

	// Procedures:
// ---------------------------------------------------------------------------------------------------------------------------------------

	// Ray memory manipulation:
	int allocate_rays ();
	int deallocate_rays ();
	
	// Upwind (and eventually, downwind) interpolation:
	int compute_upwind_points(int g, int m, Point ** Grid);
	int interpolate_simple (int g, int m, Point ** Grid);
	
	// Allocation/deallocation of the transport coefficients:
	int allocate_pqr ();
	int deallocate_pqr ();
	int allocate_rxry ();
	int deallocate_rxry ();
	int allocate_rs ();
	int deallocate_rs ();
	
	int allocate_r_full ();
	int deallocate_r_full ();
	
	// Computation of the transport coefficients:
	int compute_pqr (int g, int m); // This is the one which works best. See the notes for the specifications
	int compute_pqr_4th_direction (int g, int m); // The same function as above except contributions are asymmetric
	
	int compute_pqr_explicit_derivatives (int g, int m); // This is the one which works good. See the notes for the specifications
	
	int compute_pqr_explicit_full (int g, int m); // This is the one which probably works the best.
	
	// Intensity allocation/deallocation and manipulation:
	int allocate_I ();
	int deallocate_I ();
	int formal_solution_variant4 (int g, int m, Point ** Grid);
	
	// Formal solution related stuff:
	int formal_solution (int g, int m);
	int formal_solution_explicit_derivatives (int g, int m, Point ** Grid, int iter);
	int formal_solution_explicit_derivatives_factors (int g, int m, Point ** Grid, int iter);
	int formal_solution_explicit_full (int g, int m, Point ** Grid, int iter);
	int formal_solution_explicit_full_factors (int g, int m, Point ** Grid, int iter);
	
	int update_intensity (int g, int m);
	
	// Derivatives computation: NOT NEEDED FOR MOST OF THE APPROACHES! 
	int compute_source_derivatives (Point ** Grid);
	int compute_source_derivatives_mod (Point ** Grid);
	int modify_derivatives_asymmetric (Point ** Grid);
	
	// Source function computation:
	
	int compute_S_simple ();
	int compute_S_simple_explicit_derivatives (Point ** Grid);
	int compute_S_simple_explicit_derivatives_factors (Point ** Grid);
	int compute_S_simple_explicit_full (Point ** Grid);
	int add_a_and_c(Point ** Grid);
	int modify_a (Point ** Grid);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////// EXPERIMENTAL FUNCTIONS ADDED AFTER 16.09. 2013. ////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Variant 1: qS_u is included in the iteration factor in directions one and two, from k-th iteration k = 0,1,....
	
	int compute_pqr_explicit_full_v1 (int g, int m);
	int formal_solution_explicit_full_v1 (int g, int m, Point ** Grid, int iter);
	
	// Variant 2: Same as currently working FBILI but with iteration factor in three directions. Idea is to use it after few FBILI iterations. 
	
	//int compute_pqr_explicit_full_v2 (int g, int m); - Not needed, as the computation of the coefficients is the same...
	int formal_solution_explicit_full_v2 (int g, int m, Point ** Grid, int iter); // This needs some slight modifications.
	
	// Another take on the classic variant
	
	int formal_solution_classic (int g, int m, Point ** Grid, int iter);
	int compute_pqr_classic (int g, int m);
	int compute_S_classic(Point ** Grid);

	// Some functions needed for periodic boundaries:

	int compute_periodic_boundaries (int g, int m);
		
// -------------------------------------------------------------------------------------------------------------------
};

#endif

