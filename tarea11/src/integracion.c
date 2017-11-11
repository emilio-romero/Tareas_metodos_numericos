#include "integracion.h"
//funcion auxiliar

double Romberg(double a, double b, int n){
  double h=(b-a)/pow(2,(double)n);

return 0.0;}

double recursive_Trapezio(double(*f)(double),int n, double a, double b){
  double ranterior=0.0; 
  double **rnn=(double**)malloc((n+1)*sizeof(double*));
  for(int k=0;k<=n;k++) rnn[k]=(double*)malloc((n+1)*sizeof(double));
  
  for(int i=0;i<=n;i++){
    rnn[i][0]=rt_Ri(f,i,ranterior,a,b);
    ranterior=rnn[i][0];
  }
  for(int j=1;j<=n;j++){
    for(int i=j;i<=n;i++){
      rnn[i][j]=extrapolation_Richardson(j,rnn[i][j-1],rnn[i-1][j-1]);;
    }
  }
  printf("%g \n",rnn[n][n]);
  for(int k=0;k<=n;k++) free(rnn[k]);
  free(rnn);

return 0.0;}

double rt_Ri(double(*f)(double), int ii, double rim1,double a, double b){
 if(ii>0){
    double h=(b-a)/pow(2,(double)ii);
    double st=0;
    int ls=(int)pow(2,(double)ii-1.0);
    for(int i=0;i<ls;i++){
      st+=f(a+(2.0*((double)i+1.0)-1.0)*h);
    }
    st=st*h + (rim1/2.0);
    return(st);
  } else{
    double h=(b-a)/2.0;
    return(h*(f(a)+f(b)));
  }
}

double extrapolation_Richardson(int m, double rn, double rnm1){
  double cociente=1.0/(pow(4.0,(double)m)-1.0);
  return(rn+(rn-rnm1)*cociente);
}

double f1(double x){
  return(sin(2.0*M_PI*x));
}

double f2(double x){
  return(4.0*x*x*x-2.0*x+1.0);
}
