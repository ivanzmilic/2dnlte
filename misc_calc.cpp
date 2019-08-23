#include "all.h"

using namespace std;

double voigt(double a, double v){

	int n;
	
	alglib::real_1d_array xa;
	alglib::real_1d_array wa;
	alglib::ae_int_t info;
	alglib::ae_int_t N=10;

	// Sad izracunamo N, u zavisnosti od duzine:

	double NN = aps(v)/0.25;
	double intpart;
	modf(NN, &intpart);
	N = intpart;

	alglib::gqgenerategausslegendre(N, info, xa, wa);

	xa.setlength(N);
	wa.setlength(N);

	v = aps(v);

	double * x = xa.getcontent();
	double * w = wa.getcontent();

	for (n=0; n<N; n++){
		x[n] = x[n] * v/2 + v/2;
		w[n] = w[n] * v/2;
	}

	double I = 0.0;
	for (n=0; n<N; n++)
		I+= pow(e, x[n]*x[n] - v*v) * sin(2*a*(v-x[n])) * w[n];

	return (pow(e, a*a - v*v) * alglib::errorfunctionc(a) * cos (2*a*v) + 2.0 * pisqrt * I) * 1/pisqrt;
}

double interpolate (double x1, double x2, double x3, double y1, double y2, double y3, double x){
	

	double result;
	double w1, w2, w3;

	if (x1 == x2){
		w2 = (x3 - x) / (x3 - x2);
		w3 = (x2 - x) / (x2 - x3);
		return w2 * y2 + w3 * y3;
	}

	else if (x>x2) {

		if (aps(x-x2) / (x3-x2) < 1e-2) return y2;

		w1 = (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1);
		w2 = (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2);
		w3 = (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3);

		result = w1 * y1 + w2 * y2 + w3 * y3;

		// This part regulates what happens if there is an unrealistic extremum. If there is we return linear interpolation in between last two:
		if ((result - y2) * (result - y3) >= 0){
			w2 = (x3 - x) / (x3 - x2);
			w3 = (x2 - x) / (x2 - x3);
			return w2 * y2 + w3 * y3;
		}

		else return result;
	}

	else {

		if (aps(x-x2) / (x2-x1) < 1e-2) return y2;

		w1 = (x2 - x) * (x3 - x) / (x2 - x1) / (x3 - x1);
		w2 = (x1 - x) * (x3 - x) / (x1 - x2) / (x3 - x2);
		w3 = (x1 - x) * (x2 - x) / (x1 - x3) / (x2 - x3);

		result = w1 * y1 + w2 * y2 + w3 * y3;

		// This part regulates what happens if there is an unrealistic extremum. If there is we return linear interpolation in between last two:
		if ((result - y1) * (result - y2) >= 0){
			w1 = (x2 - x) / (x2 - x1);
			w2 = (x1 - x) / (x1 - x2);
			return w1 * y1 + w2 * y2;
		}

		else return result;
	}
}

double interpolate_bezier (double x1, double x2, double x3, double y1, double y2, double y3, double x){
	

	double result;
	double w1, w2, w3;

	if (x1 == x2){ // Then linear
		w2 = (x3 - x) / (x3 - x2);
		w3 = (x2 - x) / (x2 - x3);
		return w2 * y2 + w3 * y3;
	}

	if (x> x2){
		
		double h1 = x2 - x1;
		double h2 = x3 - x2;
		double d1 = (y2 - y1) / h1;
		double d2 = (y3 - y2) / h2;
		double a = 0.3333 * (1.0 + h2 / (h1 + h2));
		double yp;
		double C;
		double u = (x - x2) / (x3 - x2);
		
		if (d1 * d2 <= 0)
			yp = 0.0;
		else 
			yp = d1 * d2 / (a * d2 + (1.0 - a) * d1);
		
		C = y2 + h2 * 0.5 * yp;
		
		return (1.0 - u) * (1.0 - u) * y2 + u * u * y3 + 2.0 * u * (1.0 - u) * C;  

	}
	
	double h1 = x2 - x1;
	double h2 = x3 - x2;
	double d1 = (y2 - y1) / h1;
	double d2 = (y3 - y2) / h2;
	double a = 0.3333 * (1.0 + h2 / (h1 + h2));
	double yp;
	double C;
	double u = (x - x1) / (x2 - x1);
	
	if (d1 * d2 <= 0)
		yp = 0.0;
	else 
		yp = d1 * d2 / (a * d2 + (1.0 - a) * d1);
	
	C = y2 - h1 * 0.5 * yp;
	
	return (1.0 - u) * (1.0 - u) * y1 + u * u * y2 + 2.0 * u * (1.0 - u) * C;  
	
}


double aps (double x){
	
	if (x < 0) return -x;
	else return x;
}


//double w1 = (x - x2) * (x - x3) / (x1 - x2) / (x1 - x3);
		//double w2 = (x - x1) * (x - x3) / (x2 - x1) / (x2 - x3);
		//double w3 = (x - x1) * (x - x2) / (x3 - x2) / (x3 - x1);
		
		//return w1 * y1 + w2 * y2 + w3 * y3;
