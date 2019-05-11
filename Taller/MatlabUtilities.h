
#ifndef MATLABUTILITIES_H_INCLUDED
#define MATLABUTILITIES_H_INCLUDED
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include "rpoly.h"
#define NELEMS(x) (sizeof(x)/sizeof((x)[0]))
#include <stdlib.h>

long double modulo(long double a, long double b);

 long double Norma( long double v[]);

 long double dot(int dim,  long double v1[dim],  long double v2[dim]);

int Sing( long double x);

void cross( long double v1[],  long double v2[],  long double cr[3]);

 long double det(int dimension,  long double matrix[dimension][dimension]);

void multiplicacion(int r1, int c1, int r2, int c2,  long double a[r1][c1],  long double b[r2][c2],  long double c[r1][c2]);

void transpuesta(int filas, int columnas,  long double a[filas][columnas],  long double transpuesta[columnas][filas]);

void zeros(int m, int n,  long double matriz[m][n]);

bool isReal( long double parteImaginaria);

int raicesPolinomiales( int degree,  long double op[], long double zeror[degree],  long double zeroi[degree], long double resultado[20]);

 long double fix( long double x);

bool fequal( long double a,  long double b);
#endif // MATLABUTILITIES_H_INCLUDED
