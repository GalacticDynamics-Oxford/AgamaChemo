#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string.h>
int *imatrix(int);
int **imatrix(int,int);
int ***imatrix(int,int,int);
void delmatrix(int**,int);
void delmatrix(int***,int,int);
double *dmatrix(int);
double **dmatrix(int,int);
double ***dmatrix(int,int,int);
double ****dmatrix(int,int,int,int);
void delmatrix(double**,int);
void delmatrix(double***,int,int);
void delmatrix(double****,int,int,int);
float *fmatrix(int);
float **fmatrix(int,int);
float ***fmatrix(int,int,int);
float ****fmatrix(int,int,int,int);
float *****fmatrix(int,int,int,int,int);
void delmatrix(float**,int);
void delmatrix(float***,int,int);
void delmatrix(float****,int,int,int);
void delmatrix(float*****,int,int,int,int);

void compress(FILE*,float* ,int);
void compress(FILE*,double*,int);
void compress(FILE*,float** ,int,int);
void compress(FILE*,double**,int,int);
void compress(FILE*,float*** ,int,int,int);
void compress(FILE*,double***,int,int,int);
void compress(FILE*,float**** ,int,int,int,int);
void compress(FILE*,double****,int,int,int,int);
void compress(FILE*,std::vector<float >&);
void compress(FILE*,std::vector<double>&);

bool get(FILE*,float* ,int);
bool get(FILE*,double*,int);
bool get(FILE*,float** ,int,int);
bool get(FILE*,double**,int,int);
bool get(FILE*,float*** ,int,int,int);
bool get(FILE*,double***,int,int,int);
bool get(FILE*,float**** ,int,int,int,int);
bool get(FILE*,double****,int,int,int,int);
bool get(FILE*,float***** ,int,int,int,int);
bool get(FILE*,double*****,int,int,int,int);

#include "grid2.h"
#include "grid3.h"
#include "grid4.h"