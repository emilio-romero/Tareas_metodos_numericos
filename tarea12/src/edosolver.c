#include "edosolver.h"

void edo_trapezio(int n, double li, double ls, double yini, double dyini){
double **mati=(double**)malloc(2*sizeof(double*));
double *vecd=(double*)malloc(2*sizeof(double));
double *ysol=(double*)malloc(2*sizeof(double));
double *yip1=(double*)malloc(2*sizeof(double));
FILE *salida; 
salida=fopen("Solucion.dat","w");
for(int i=0;i<2;i++) mati[i]=(double*)malloc(2*sizeof(double));
  double h=(ls-li)/(double)n;
  double xi,xip1,error=0;
    mati[0][0]=1;
    mati[0][1]=-h/2.0;
    mati[1][0]=h/2.0;
    ysol[0]=yini;
    ysol[1]=dyini;
    fprintf(salida,"%lf ",li);
    fprintf(salida,"%lf\n",ysol[0]);
  for(int i=1;i<n;i++){
    xi=li+h*(double)(i-1);
    xip1=li+h*(double)(i);
    mati[1][1]=1.0-h*xip1/2.0;
    vecd[0]=h*ysol[1]/2.0+ysol[0];
    vecd[1]=ysol[1]+h*(-ysol[0]+xi*ysol[1]-2.0*(xi*cos(xi)+xip1*cos(xip1)))/2.0;
    /*Obtencion de Y_{i+1}*/
    yip1=factLU(mati,vecd,2); 
    ysol[0]=yip1[0];
    ysol[1]=yip1[1];
    fprintf(salida,"%lf ",xip1);
    fprintf(salida,"%lf\n",ysol[0]);

    if(fabs(ysol[0]-(xip1+2.0*sin(xip1)))>error)
      error=fabs(ysol[0]-(xip1+2.0*sin(xip1)));
  }
  printf("Ultima solucion: %g \n",ysol[0]);
  printf("Error: %g\n",error);
  fclose(salida);
for(int i=0;i<2;i++) free(mati[i]);
free(mati);
free(vecd);
free(ysol);
}
