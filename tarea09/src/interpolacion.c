#include "interpolacion.h"

void Vandermonde(double **A, double *x, int n){
  for(int i=0;i<n;i++){
    A[i][0]=1.0;
    A[i][1]=x[i];
    for(int j=2;j<n;j++){
      A[i][j]=pow(x[i],j);
    }
  }
}
double *polinomial(double *x0, double *y0, int n){
double *c=(double*)malloc(n*sizeof(double));
double **A=(double**)malloc(n*sizeof(double*));
for(int i=0;i<n;i++) A[i]=(double*)malloc(n*sizeof(double));

c=factLU(A,c,n);

return(c);}
