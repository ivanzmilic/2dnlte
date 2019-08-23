//#ifndef GLOBALS_H
//#define GLOBALS_H

class Globals{
private:
	static double * x; // This is frequency grid in Doppler width units (for now)
	static int N; // Resolution of the frequency grid

	static int NX; // Geometrical resolution in X
	static int NY; // Geometrical resolution in Y
	
// ------------------------------------------------------------------------------------

	// Cylindric angles. In case we are not using different grid for each point:

	static int NP;
	static int NT;
	static int * NP_Carlson; // This follows from the paper of Carlson, as the number of angles along the azimuth DEPENDS on the particular latitude we are talking about.
	static double * phi;
	static double ** phi_C; // As in "phi_Carlson", just shorter
	static double ** w_C; //  For "weights_Carlson"
	static double * theta;
	static double * sinphi;
	static double * cosphi;
	static double * sintheta;
	static double * costheta;
	static double * wphi;
	static double * wtheta;
	
	// Now there are sines and cosines for carlson angle set:
	
	static double ** sinphi_C;
	static double ** cosphi_C;
	
	// Now there is an additional beauty...
	// These are not REAL (i.e. ``problem'' angles as Van Noort et al would say). These are local angles, that is theta and phi are actually gamma and varphi
	// or whatever. Now we will define here a Jacobian which will compute what needs to be computed in order to swap integration:
	
	static int type; // Type of the profile: 1 - Doppler; 2 - Voigt
	static double * profile;
	static double * wx;

	static double * X;
	static double * Y;
	
	static double **** Psi; // This is Psi matrix from Anusha & Nagendra (2011b)
	
	// Additional useful quantities:
	
	static double w; // Over-relaxation factor for SOR/SSOR techniques
	static int factor; // A quantity which determines if the iteration factor is used and how.

public:
	static void allocate_x (int N_in);
	static void compute_x ();
	static double getx(int n);
	static int getN();
	static void compute_profile();
	static double getprofile(int n);
	static double getwx(int n);

// Cylindrical angles

	static void setNP (int NP_in);
	static void setNT (int NT_in);
	static int getNP ();
	static int getNT ();
	static void compute_theta ();
	static void compute_phi ();
	
	static int compute_angles_Carlson ();
	static int compute_angles_Carlson_10 ();
	static int compute_angles_simple (int NTnew, int NPnew);
	static int compute_angles_pseudogauss (int NTnew, int NPnew);
	static int compute_angles_1 ();
	
	static double gettheta(int g);
	static double getsintheta (int g);
	static double getcostheta (int g);
	
	// Now, get functions for phi, cosphi, and sinphi in Carlson set:
// ------------------------------------------------------------------
	static double getphi_C (int g, int m);
	static double getcosphi_C (int g, int m);
	static double getsinphi_C (int g, int m);
	static double getw_C (int g, int m);
	static int getNP_C (int g);
// ------------------------------------------------------------------

	static void setNX (int NXin);
	static int getNX();

	static void setNY (int NYin);
	static int getNY ();

	static void deallocate_all();

	// Some geometry related things:
	static int set_geometry (Point ** Grid);
	static double * get_X_axis ();
	static double * get_Y_axis ();
	static double getX (int i);
	static double getY (int j);
	
	// And real matrix for working ni reduced bases
	
	static void compute_Psi ();
	static double getPsi (int l, int lp, int g, int m);
	
	// Manipulation functions for these useful additional quantities:
	
	static void set_w (double intw);
	static double get_w ();
	static void set_factor (int infactor);
	static int get_factor();

};

//#endif
