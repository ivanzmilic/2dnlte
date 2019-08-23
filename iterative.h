#ifndef ITERATIVE_H
#define ITERATIVE_H

double FBILI_2by2 (Point ** Grid, int N_iterations, double delta);

double FBILI_2by2_explicit_derivatives (Point ** Grid, int N_iterations, double delta);

double SSOR_explicit_derivatives (Point ** Grid, int N_iterations, double delta);

double Jacobi_explicit_derivatives (Point ** Grid, int N_iterations, double delta);

double GS_explicit_derivatives (Point ** Grid, int N_iterations, double delta);

double FBILI_explicit_derivatives (Point ** Grid, int N_iterations, double delta);

// ------------------------------------------------------------------------------------------------------------------------------------------------------

// Now we shall put here FULLY explicit versions, i.e. versions which explicitly take contributions of the S_u to the intensity as well as the derivatives.
// Contrary to previous coding we shall keep old names for r1-r4, and add r5-r8 as the "corner" source functions.

double Jacobi_explicit_full (Point ** Grid, int N_iterations, double delta);

double GS_explicit_full (Point ** Grid, int N_iterations, double delta);

double SSOR_explicit_full (Point ** Grid, int N_iterations, double delta);

double FBILI_explicit_full_2factors (Point ** Grid, int N_iterations, double delta);

double FBILI_explicit_full_2by2 (Point ** Grid, int N_iterations, double delta);

// -------------------------------------------------------------------------------------------------------------------------------------------------------

// Now we will try the old FBILI procedure:

double FBILI_upwind_derivative (Point ** Grid, int N_iterations, double delta);

// Read S from file, and other "help" functions:

int read_S_from_file (Point ** Grid);

// Re-attempt @ elmination scheme

double Jacobi_classic (Point ** Grid, int N_interations, double delta);



#endif

