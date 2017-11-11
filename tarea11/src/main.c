#include <stdio.h> 
#include "integracion.h"
double mf(double x){
  return x;
}
int main(){

recursive_Trapezio(mf,6,0,10);
recursive_Trapezio(f1,6,0,M_PI);
recursive_Trapezio(f2,6,0,M_PI);

printf("Su programa ha terminado\n");
return 0;}
