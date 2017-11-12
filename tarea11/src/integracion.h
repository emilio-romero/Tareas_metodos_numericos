#ifndef INTEGRACION_H 
#define INTEGRACION_H 
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 
double rt_Ri(double(*f)(double),int ii, double rim1, double a, double b);
double Romberg(double(*f)(double),int n,double a,double b,double(*F)(double));
double extrapolation_Richardson(int m, double rn, double rnm1);
double f1(double x);
double f2(double x);
double F1(double x);
double F2(double x);
#endif 
