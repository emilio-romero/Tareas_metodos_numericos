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
int writeData(double **A, int nr, int nc, char *cfile);
/*Tarea 10*/
double Ni(double xk, double h, double zk);
double Nip1(double xk, double h, double zk);
double evalNi(double **data, double x, int ii, double h);
double evalNip1(double **data, double x, int ii, double h);
double coefa(double **data, int m, double h, double lambda, int k);
double coefb(double **data, int m, double h, double lambda, int k);
double coefc(double **data, int m, double h, double lambda, int k);
double coefs(double **data, int m, double h,  int k);
double coeft(double **data, int m, double h,  int k);
double funPhi(double *phis, double x, int n, double **data, double h,int ii);
void diffMatrix(double **A,double **data,double *b,int m,int n,double h, double lambda);
void elementosFinitos(char *cfile, int n, double lambda, char *outfile);

void splineCubico(char *cfile, int m,char *outfile);
void scMatrix(double **td, double *d, double **data, int n);
void coefM(double *M, double **data, int n);
double evalSC(double *M, double **data, double x, int n);
#endif

