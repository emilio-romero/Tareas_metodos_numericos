#include <stdio.h> 
#include "integracion.h"
double mf(double x){
  return x;
}
int main(){
printf("n=4\n");
printf("f(x)=sin(2pi x)\n");
Romberg(f1,4,0,M_PI,F1);
printf("f(x)=4x^3-2x+1\n");
Romberg(f2,4,0,M_PI,F2);
printf("n=6\n");
printf("f(x)=sin(2pi x)\n");
Romberg(f1,6,0,M_PI,F1);
printf("f(x)=4x^3-2x+1\n");
Romberg(f2,6,0,M_PI,F2);
printf("n=10\n");
printf("f(x)=sin(2pi x)\n");
Romberg(f1,10,0,M_PI,F1);
printf("f(x)=4x^3-2x+1\n");
Romberg(f2,10,0,M_PI,F2);

printf("Su programa ha terminado\n");
return 0;}
