#ifndef MATLABUTILITIES_H_INCLUDED
#define MATLABUTILITIES_H_INCLUDED
#include <math.h>

double Norma(double v[]);

double dot(int dim, double v1[dim], double v2[dim]);

int Sing(double x);

double* cross(double v1[], double v2[]);

double det(int dimension, double matrix[dimension][dimension]);

void multiplicacion(int dimension, double a[dimension][dimension], double b[dimension][dimension], double c[dimension][dimension]);

void transpuesta(int dimension, double a[dimension][dimension],  double transpuesta[dimension][dimension]);

void zeros(int m, int n, double matriz[m][n]);
#endif // MATLABUTILITIES_H_INCLUDED
