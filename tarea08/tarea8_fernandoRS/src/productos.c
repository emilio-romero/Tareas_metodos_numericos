#include "productos.h"

double prodPunto(double *a, double *b, int n){
  double suma=0;
  for(int i=0;i<n;i++){
    suma+=a[i]*b[i];
  }

return(suma);}

void prodMatVec(double *out, double **A, double *v, int n){
  double suma; 
  for(int i=0;i<n;i++){
    suma=0;
    for(int k=0;k<n;k++){
      suma+=A[i][k]*v[k];
    }
    out[i]=suma;
  }
}

void sumVec(double *out, double *a, double *b, int n){
  for(int i=0;i<n;i++){
    out[i]=a[i]+b[i];
  }
}

void prodSclVec(double *out, double a, double *v, int n){
  for(int i=0;i<n;i++){
    out[i]=a*v[i];
  }
  }
void transLinear(double *out, double a, double *b, double *v, int n){
  for(int i=0;i<n;i++){
    out[i]=v[i]+a*b[i];
  }

}
