#ifndef INTERPOLACION_H
#define INTERPOLACION_H
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include "solucionadores.h"
#include "lectura.h"
void Vandermonde(double **A, double *x, int n);
double *polinomial(double *x0, double *y0, int n);
double evalPolinomio(double *c, double x, int n);
void interpolaPolinomio(double *x0, double *y0, double xi, double xf, int n);
double Lagrange(double *x0, double x, int n, int i);
double evalLagrangep(double *x0, double *y0, double x, int n);
void interpolaLagrange(double *x0, double *y0, double xi, double xf, int n);
double Newtonp(double *x0, double x, int n, int i);
void diferenciasDivididas(double **A, double *x0, double *y0, int n);
double evalNewtonp(double **A, double *x0, double x, int n);
void interpolaNewton(double *x0, double *y0, double xi, double xf, int n);
#endif 
