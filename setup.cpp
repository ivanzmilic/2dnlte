#include "all.h"

using namespace std;

int geometry (Point ** Grid){

	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int i,j;
	
	double T_X = 1E2;
	double T_Y = 1E5;
	
	// Enumerate points:

	for (i=0; i<NX; i++)
		for (j=0; j<NY; j++){
			Grid[i][j].i = i;
			Grid[i][j].j = j;
		}

	// There can be an option to load manually grid from a file, but we are not going to do it. We are going to make few grids, and swap manually between them. 

	int grid_no = 1;

	if (grid_no == 1){ // Log spatial grid with fine spacing toward the edges, as in Auer & Paletou (1994)
		
		double Xstep;
		double Ystep;
		
		double smallestexp = - 3;

		int symmetric_Y = 0;

		int centerX = (NX-1) / 2;
		int centerY = (NY-1) / 2;

		Xstep = (log10(T_X/2.0) - smallestexp) / (centerX-1);
		
		if (symmetric_Y == 1)
			Ystep = (log10(T_Y/2.0) - smallestexp) / (centerY-1);

		else 
			Ystep = (log10(T_Y) - smallestexp) / (NY-1);
		
		for (i=0; i<centerX+1; i++)
			for (j=0; j<NY; j++){
				Grid[i][j].X = pow(10, smallestexp + (i-1) * Xstep);
				if (i==0) 
					Grid[i][j].X = 0.0;
			}


		
		for (i=centerX+1; i<NX; i++)
			for (j=0; j<NY; j++)
				Grid[i][j].X = T_X - Grid[NX-1-i][j].X;
		

		if (symmetric_Y == 1){
			
			for (i=0; i<NX; i++)
				for (j=0; j<centerY+1; j++){
					Grid[i][j].Y = pow(10.0, smallestexp + (j-1) * Ystep);
					if (j==0) 
						Grid[i][j].Y = 0;
				}
		
			for (i=0; i<NX; i++)
				for (j=centerY+1; j<NY; j++)
					Grid[i][j].Y = T_Y - Grid[i][NY-1-j].Y;
		}

		else {

			for (i=0; i<NX; i++)
				for (j=0; j<NY; j++){
					Grid[i][j].Y = pow(10.0, smallestexp + (j-1) * Ystep);
					if (j==0) 
						Grid[i][j].Y = 0;
				}

		}

		int equidistant_X = 1;

		if (equidistant_X == 1){

			Xstep = T_X / (NX-1);

			for (i=0; i<NX; i++)
				for (j=0; j<NY;j++)

					Grid[i][j].X = i *  Xstep;

		}
	}

	// Now copy relevant stuff into the globals:
	
	Globals::set_geometry(Grid);
// -----------------------------------------------------------------------------------------------------------------------------

	return 0;
}

int RT_parameters (Point ** Grid){

	int NX = Globals::getNX();
	int NY = Globals::getNY();
	int i,j;

	for (i=0; i<NX; i++)
		for (j=0; j<NY; j++){
			Grid[i][j].B = 1.0;// + 0.0 * cos(Grid[i][j].X * 2 * pi / 1.*Globals::getX(NX-1));
			Grid[i][j].Chi = 1.0; // We will write the procedure for interpolation regardless of whatever happens.
			Grid[i][j].ep = 1E-4;
			Grid[i][j].S = 1; //Grid[i][j].B; // Initial vallue
			Grid[i][j].dS = 0;
			Grid[i][j].J = 0;
			Grid[i][j].alo = 0;
			Grid[i][j].BB = 0;
			Grid[i][j].C1 = 0;
			Grid[i][j].C2 = 0;
			Grid[i][j].C3 = 0;
			Grid[i][j].C4 = 0;
			Grid[i][j].AA = 0;
			Grid[i][j].CX = 0;
			Grid[i][j].CY = 0;
			
			// In case you want to explicitly treat also upwind source function:
			
			Grid[i][j].C5 = 0;
			Grid[i][j].C6 = 0;
			Grid[i][j].C7 = 0;
			Grid[i][j].C8 = 0;
			
			// In case you separately treat AA coefficients:
			
			Grid[i][j].AA1 = 0;
			Grid[i][j].AA2 = 0;
			Grid[i][j].AA3 = 0;
			Grid[i][j].AA4 = 0;
			
			Grid[i][j].BB1 = 0;
			Grid[i][j].BB2 = 0;
			Grid[i][j].BB3 = 0;
			Grid[i][j].BB4 = 0;
		}
	//Grid[50][40].B+=0.05;

	return 0;
}

