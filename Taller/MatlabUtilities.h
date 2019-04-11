#ifndef MATLABUTILITIES_H_INCLUDED
#define MATLABUTILITIES_H_INCLUDED


double Norma(double v[]);

double dot(double v1[], double v2[]);

int Sing(double x);

double* cross(double v1[], double v2[]);

double det(int dimension, double matrix[dimension][dimension]);

void multiplicacion(double a[3][3], double b[3][3], double c[3][3]);

void transpuesta(double a[3][3],  double transpuesta[3][3]);
#endif // MATLABUTILITIES_H_INCLUDED
