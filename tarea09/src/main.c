#include <stdio.h> 
#include <stdlib.h> 
#include "interpolacion.h"
int main(int argc, char *argv[]){
char archivo[30]; 
readParams(argc, argv, archivo);
int nr, nc; 
double **mydata; 
double *xs, *ys;
mydata=readMatrix(archivo,&nr,&nc);
printf("row: %d , cols: %d\n",nr,nc);
xs=(double*)malloc(nr*sizeof(double));
ys=(double*)malloc(nr*sizeof(double));
for(int i=0;i<nr;i++){
  xs[i]=mydata[i][0];
  ys[i]=mydata[i][1];
}
interpolaPolinomio(xs,ys,mydata[0][0],mydata[nr-1][0],nr);
interpolaLagrange(xs,ys,mydata[0][0],mydata[nr-1][0],nr);
interpolaNewton(xs,ys,mydata[0][0],mydata[nr-1][0],nr);

freeMatrix(mydata);
printf("Su programa ha terminado\n");
return 0;}
