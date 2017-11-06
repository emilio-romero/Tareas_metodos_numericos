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
/*Metodo de los elementos finitos para interpolacion*/

double Ni(double x, double h, double zk){
  double aux; 
  aux=1.0-(1.0/h)*(x-zk);
return aux;}

double Nip1(double x, double h, double zk){
  double aux; 
  aux=(x-zk)/h;
return aux;}

double evalNi(double **data,double x,int ii, double h){
  double ls, li, res; 
  ls=data[0][0]+h*(double)(ii+1);
  li=data[0][0]+h*(double)ii;
    if(x>=li && x<ls){
      res=Ni(x,h,li);
    } else{
    res=0.0;
    }
return res;}

double evalNip1(double **data, double x, int ii, double h){
  double ls, li, res; 
  ls=data[0][0]+h*(double)(ii);
  li=data[0][0]+h*(double)(ii-1);
    if(x>=li && x<ls){
      res=Nip1(x,h,li);
    } else{
    res=0.0;
    }

return res;}

double coefa(double **data, int m, double h, double lambda, int k){
  double suma=0.0;
  for(int i=0;i<m;i++){
    suma+=evalNi(data,data[i][0],k,h)*evalNi(data,data[i][0],k,h);
  }
  suma+=lambda/h;

return suma;}
double coefb(double **data, int m, double h, double lambda, int k){
  double suma=0.0;
  for(int i=0;i<m;i++){
    suma+=evalNip1(data,data[i][0],k+1,h)*evalNi(data,data[i][0],k,h);
  }
  suma=suma-lambda/h;
return suma;}
double coefc(double **data, int m, double h, double lambda, int k){
  double suma=0.0;
  for(int i=0;i<m;i++){
    suma+=evalNip1(data,data[i][0],k+1,h)*evalNip1(data,data[i][0],k+1,h);
  }
  suma+=lambda/h;
return suma;}
double coefs(double **data, int m, double h, int k){
  double suma=0.0;
  for(int i=0;i<m;i++){
    suma+=data[i][1]*evalNi(data,data[i][0],k,h);
  }
return suma;}
double coeft(double **data, int m, double h, int k){
  double suma=0.0;
  for(int i=0;i<m;i++){
    suma+=data[i][1]*evalNip1(data,data[i][0],k+1,h);
  }
return suma;}

void diffMatrix(double **A,double **data,double *b,int m,int n,double h,double lambda){
  for(int j=1;j<n-1;j++){
    A[0][j]=coefb(data,m,h,lambda,j); //derecha  
    A[1][j]=coefc(data,m,h,lambda,j)+coefa(data,m,h,lambda,j);//centro
    A[2][j]=coefb(data,m,h,lambda,j);//izquierda
    b[j]=coeft(data,m,h,j-1)+coefs(data,m,h,j);
  }
  A[0][0]=coefb(data,m,h,lambda,0);
  A[1][0]=coefa(data,m,h,lambda,0);
  A[2][0]=0.0;
  A[0][n-1]=0.0;
  A[1][n-1]=coefc(data,m,h,lambda,n-1);
  A[2][n-1]=coefb(data,m,h,lambda,n-1);
  b[0]=coefs(data,m,h,0);
  b[n-1]=coeft(data,m,h,n-1);
  //printVector(A[0],n);
  //printMatrix(A,3,n);
}

double funPhi(double *phis, double x,int n,double **data, double h, int ii){
  double aux; 
    aux=phis[ii]*evalNi(data,x,ii,h)+phis[ii+1]*evalNip1(data,x,ii+1,h);
return aux;}

void elementosFinitos(char *cfile, int n,double lambda, char *outfile){
double **mydata; 
int nr, nc, pos; 
  mydata=readMatrix(cfile,&nr,&nc);
  double h=(mydata[nr-1][0]-mydata[0][0])/((double)n);
  double xi,yi,zk;
  double lims, limi;
  double **td=(double**)malloc(3*sizeof(double*));
  double *b=(double*)malloc(n*sizeof(double));
  double *phi=(double*)malloc(n*sizeof(double));
  for(int k=0;k<3;k++) td[k]=malloc(n*sizeof(double));
  diffMatrix(td,mydata,b,nr,n,h,lambda);
  phi=tridiag(td[2],td[1],td[0],b,n); 
  //printVector(phi,n);
 double **interData=(double**)malloc((nr-1)*sizeof(double*));
  for(int i=0;i<nr-1;i++) interData[i]=(double*)malloc(2*sizeof(double));
  for(int i=0;i<nr-1;i++){
    interData[i][0]=mydata[i][0];
    for(int j=0;j<n;j++){
      limi=mydata[0][0]+h*(double)j;
      lims=mydata[0][0]+h*(double)(j+1);
      if(mydata[i][0]>=limi && mydata[i][0]<lims){
        pos=j;
        break;
      }
    }
    interData[i][1]=funPhi(phi,mydata[i][0],n,mydata,h,pos); 
  }
  writeData(interData,nr-1,2,outfile);
for(int j=0;j<nr-1;j++) free(interData[j]);
free(interData);
 
  for(int k=0;k<3;k++) free(td[k]);
  free(td);
  free(b);
  free(phi);
  freeMatrix(mydata);
}
/*Interpolacion por splines cubicos*/

void scMatrix(double **td, double *d, double **data, int n){
  double  hi, hip1;
  int id;
  for(int i=0;i<n-2;i++){
    id=i+1;
    hi=data[id][0]-data[id-1][0];
    hip1=data[id+1][0]-data[id][0];
    td[0][i]=hip1/(hi+hip1);//derecho
    td[1][i]=2; //central
    td[2][i]=hi/(hi+hip1);//izquierdo
    d[i]=6/(hi+hip1)*((data[id+1][1]-data[id][1])/hip1-(data[id][1]-data[id-1][1])/hi);
  }
  td[0][n-3]=0;
  td[1][0]=2;
  td[2][0]=0;
}

void coefM(double *M, double **data,int n){
  M[0]=0.0; M[n-1]=0.0;
  double *mr=(double*)malloc((n-2)*sizeof(double));
  double *b=(double*)malloc((n-2)*sizeof(double));
  double **td=(double**)malloc(3*sizeof(double*));
  for(int i=0;i<3;i++) td[i]=(double*)malloc((n-2)*sizeof(double));
  scMatrix(td,b,data,n);
  mr=tridiag(td[2],td[1],td[0],b,n-2);
  for(int i=1;i<n-1;i++){
    M[i]=mr[i-1];
  }
for(int i=0;i<3;i++) free(td[i]);
free(td);
free(mr);
free(b);
}

double evalSC(double *M, double **data, double x, int n){
double eval; 
double hi, cih, ci;  
double xi, xim1;
int ind;
  for(int i=0;i<n;i++){
    if(x<data[i][0]){
      ind=i;
      break;
    }
  }
xi=data[ind][0]-x;
xim1=x-data[ind-1][0];
hi=data[ind][0]-data[ind-1][0]; 
cih=data[ind-1][1]-M[ind-1]*hi*hi/6.0; 
ci=(data[ind][1]-data[ind-1][1])/hi -(hi/6.0)*(M[ind]-M[ind-1]);
eval=M[ind-1]*xi*xi*xi/(6.0*hi) +M[ind]*xim1*xim1*xim1/(6.0*hi)+ci*xim1+cih;

return eval;}

void splineCubico(char *cfile, int m, char *outfile){
double **mydata; 
  int nr, nc; 
  mydata=readMatrix(cfile,&nr,&nc);
  double *Ms=(double*)malloc(nr*sizeof(double));
  coefM(Ms,mydata,nr);
  double dx=(mydata[nr-1][0]-mydata[0][0])/(double)m;
  double **interData=(double**)malloc(m*sizeof(double*));
  for(int i=0;i<m;i++) interData[i]=(double*)malloc(2*sizeof(double));
  for(int i=0;i<m;i++){
    interData[i][0]=mydata[0][0]+dx*(double)i;
    interData[i][1]=evalSC(Ms,mydata,interData[i][0],nr); 
  }
  writeData(interData,m,2,outfile);
for(int j=0;j<m;j++) free(interData[j]);
free(interData);
free(Ms);
freeMatrix(mydata);
}


