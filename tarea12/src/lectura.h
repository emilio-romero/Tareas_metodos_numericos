#ifndef LECTURA_H
#define LECTURA_H 
void readParams(int argc, char *argv[],char *cfile); 
double *readVector(char *cfile, int *nr);
int writeVector(double *vec, int dim, char *cfile);
void printVector(double *vec, int dim); 
double **createMatrix(int nr, int nc); 
double **readMatrix(char *cfile, int *nr, int *nc); 
int writeMatrix(double **mat, int nr, int nc, char *cfile); 
void printMatrix(double **mat, int nr, int nc);
void freeMatrix(double **mat);
#endif 
