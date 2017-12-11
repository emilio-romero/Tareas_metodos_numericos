#include <stdio.h>
#include <stdlib.h> 
#include "raices.h" 
#include "solucionadores.h"
#include "lectura.h"
#include <float.h>
int main(int argc, char *argv[]){
char fmat[30];
char fvec[30];
readParams(argc, argv, fmat, fvec);
double **mat, *vec; 
int nr,nc;
double eps=sqrt(DBL_EPSILON);
vec=readVector(fvec,&nr);
mat=readMatrix(fmat,&nr,&nc);
double *x=(double*)malloc(nr*sizeof(double));
x=conjugateGradient(mat,vec,nr,eps);

printVector(x,nr);
free(vec); free(x);
freeMatrix(mat);
printf("\nSu programa ha terminado\n");
return 0;}
