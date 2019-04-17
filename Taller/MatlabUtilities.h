#ifndef MATLABUTILITIES_H_INCLUDED
#define MATLABUTILITIES_H_INCLUDED
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include "rpoly.h"
#define NELEMS(x) (sizeof(x)/sizeof((x)[0]))
#include <stdlib.h>

double Norma(double v[]);

double dot(int dim, double v1[dim], double v2[dim]);

int Sing(double x);

double* cross(double v1[], double v2[]);

double det(int dimension, double matrix[dimension][dimension]);

void multiplicacion(int r1, int c1, int r2, int c2, double a[r1][c1], double b[r2][c2], double c[r1][c2]);

void transpuesta(int filas, int columnas, double a[filas][columnas], double transpuesta[columnas][filas]);

double* zeros(int m, int n);

bool isReal(double parteImaginaria);

double* roots(double *op, int degree, double *zeror, double *zeroi);

double fix(double x);

#endif // MATLABUTILITIES_H_INCLUDED
