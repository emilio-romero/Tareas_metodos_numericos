#include "interpolacion.h"

void Vandermonde(double **A, double *x, int n){
  for(int i=0;i<n;i++){
    A[i][0]=1.0;
    A[i][1]=x[i];
    for(int j=2;j<n;j++){
      A[i][j]=pow(x[i],(double)j);
    }
  }
}
double *polinomial(double *x0, double *y0, int n){
double *c=(double*)malloc(n*sizeof(double));
double **A=(double**)malloc(n*sizeof(double*));
for(int i=0;i<n;i++) A[i]=(double*)malloc(n*sizeof(double));
  Vandermonde(A,x0,n);
  c=factLU(A,y0,n);
/*for(int i=0;i<n;i++){
  for(int j=0;j<n;j++){
    printf("%lf ",A[i][j]);
  } printf("\n");
}*/

for(int i=0;i<n;i++) {free(A[i]);}
free(A);

return(c);}

double evalPolinomio(double *c, double x, int n){
double px=0.0; 
  for(int i=0;i<n;i++)
    px+=c[i]*pow(x,(double)i);
return(px);}

void interpolaPolinomio(double *x0, double *y0,double xi, double xf,int n){
double *c;//=(double*)malloc(n*sizeof(double));
  c=polinomial(x0,y0,n);
double dx=(xf-xi)/(double)(4*n);
double **ndata=(double**)malloc(4*n*sizeof(double*));
  for(int i=0;i<4*n;i++){
    ndata[i]=(double*)malloc(2*sizeof(double));
    ndata[i][0]=xi+i*dx; 
    ndata[i][1]=evalPolinomio(c,ndata[i][0],n);
  }
  writeData(ndata,4*n,2,"polinomial_data.dat");
for(int i=0;i<4*n;i++) free(ndata[i]);
free(ndata);
free(c);
}

/*Polinomios de Lagrange*/
double Lagrange(double *x0, double x, int n, int i){
  double Lx=1; 
  for(int j=0;j<n;j++){
    if(j!=i)
      Lx*=(x-x0[j])/(x0[i]-x0[j]);
  }
return Lx;}

double evalLagrangep(double *x0, double *y0, double x, int n){
double px=0.0;
  for(int i=0;i<n;i++){
    px+=y0[i]*Lagrange(x0,x,n,i);
  }
return(px);}

void interpolaLagrange(double *x0, double *y0,double xi, double xf,int n){
double dx=(xf-xi)/(double)(4*n); printf("%lf\n",dx);
double **ndata=(double**)malloc(4*n*sizeof(double*));
  for(int i=0;i<4*n;i++){
    ndata[i]=(double*)malloc(2*sizeof(double));
    ndata[i][0]=xi+i*dx; 
    ndata[i][1]=evalLagrangep(x0,y0,ndata[i][0],n);
  }
  writeData(ndata,4*n,2,"poliLagrange_data.dat");
for(int i=0;i<4*n;i++) free(ndata[i]);
free(ndata);
}
/*Polinomios de Newton*/
double Newtonp(double *x0, double x, int n, int i){
 double Nx=1; 
 if(i!=0){
  for(int j=0;j<i;j++){
    Nx*=(x-x0[j]);
  }
 return Nx;
 }
 else return 1.0;
}

void diferenciasDivididas(double **A, double *x0, double *y0, int n){

  for(int i=0;i<n;i++) A[i][0]=y0[i]; 
  for(int j=1;j<n;j++){
    for(int i=0;i<n-j;i++){
      A[i][j]=(A[i+1][j-1]-A[i][j-1])/(x0[i+j]-x0[i]);
    }
  }
}

double evalNewtonp(double **A, double *x0, double x, int n){
double px=0; 
  for(int i=0;i<n;i++){
    px+=A[0][i]*Newtonp(x0,x,n,i);
  }
return(px);}

void interpolaNewton(double *x0, double *y0,double xi, double xf,int n){
double dx=(xf-xi)/(double)(4*n);
double **ndata=(double**)malloc(4*n*sizeof(double*));
double **A=(double**)malloc(n*sizeof(double*));
for(int k=0;k<n;k++) A[k]=(double*)malloc(n*sizeof(double));
diferenciasDivididas(A, x0, y0, n);
  for(int i=0;i<4*n;i++){
    ndata[i]=(double*)malloc(2*sizeof(double));
    ndata[i][0]=xi+i*dx; 
    ndata[i][1]=evalNewtonp(A,x0,ndata[i][0],n);
  }
  writeData(ndata,4*n,2,"poliNewton_data.dat");
for(int i=0;i<4*n;i++) free(ndata[i]);
free(ndata);
}

int writeData(double **A, int nr, int nc, char *cfile){
FILE *out; 
out=fopen(cfile,"w");
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      fprintf(out,"%lf ",A[i][j]);
    }
    fprintf(out,"\n");
  }
fclose(out);
return 0;}

