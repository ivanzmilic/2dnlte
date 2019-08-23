#ifndef MISC_CALC_H
#define MISC_CALC_H

double voigt(double a, double v);

double interpolate (double x1, double x2, double x3, double y1, double y2, double y3, double x);

double interpolate_bezier (double x1, double x2, double x3, double y1, double y2, double y3, double x);

double aps (double x);

#endif
