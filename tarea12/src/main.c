#include <stdio.h> 
#include "edosolver.h"
int main(){
printf("================\n");
printf("n=400\n");
printf("================\n");
edo_trapezio(400,0.0,5.0,0.0,3.0);

printf("================\n");
printf("n=4000\n");
printf("================\n");
edo_trapezio(4000,0.0,5.0,0.0,3.0);

printf("Su programa ha terminado\n");
return 0;}
