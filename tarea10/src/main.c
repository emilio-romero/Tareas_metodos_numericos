#include <stdio.h> 
#include <stdlib.h> 
#include "interpolacion.h"
int main(int argc, char *argv[]){
double lamb;
if(argc>1) lamb=atof(argv[1]);
else lamb=1.0;
elementosFinitos("datos2D_1.bin",25,lamb,"ejercicio1EF.dat");
elementosFinitos("datos2D_2.bin",25,lamb,"ejercicio2EF.dat");
splineCubico("datos2D_3.bin",100,"ejercicio1SC.dat");
splineCubico("datos2D_4.bin",100,"ejercicio2SC.dat");
printf("Su programa ha terminado\n");
return 0;}
