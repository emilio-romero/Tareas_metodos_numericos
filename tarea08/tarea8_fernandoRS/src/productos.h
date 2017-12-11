#ifndef PRODUCTOS_H
#define PRODUCTOS_H 

double prodPunto(double *a, double *b, int n);
void prodMatVec(double *out, double **A, double *v, int n);
void sumVec(double *out, double *a, double *b, int n);
void prodSclVec(double *out, double a, double *v, int n);
void transLinear(double *out, double a, double *b, double *v, int n);
#endif 
