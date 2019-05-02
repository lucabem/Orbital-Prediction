#ifndef MATLABUTILITIES_H_INCLUDED
#define MATLABUTILITIES_H_INCLUDED
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include "rpoly.h"
#define NELEMS(x) (sizeof(x)/sizeof((x)[0]))
#include <stdlib.h>

 double Norma( double v[]);

 double dot(int dim,  double v1[dim],  double v2[dim]);

int Sing( double x);

void cross( double v1[],  double v2[],  double cr[3]);

 double det(int dimension,  double matrix[dimension][dimension]);

void multiplicacion(int r1, int c1, int r2, int c2,  double a[r1][c1],  double b[r2][c2],  double c[r1][c2]);

void transpuesta(int filas, int columnas,  double a[filas][columnas],  double transpuesta[columnas][filas]);

void zeros(int m, int n,  double matriz[m][n]);

bool isReal( double parteImaginaria);

int raicesPolinomiales( int degree,  double op[], double zeror[degree],  double zeroi[degree], double resultado[20]);

 double fix( double x);

bool fequal( double a,  double b);
#endif // MATLABUTILITIES_H_INCLUDED
