#ifndef MATLABUTILITIES_H_INCLUDED
#define MATLABUTILITIES_H_INCLUDED


double Norma(double v[]);

double dot(double v1[], double v2[]);

int Sing(double x);

double* cross(double v1[], double v2[]);

double det(int dimension, double matrix[dimension][dimension]);

void multiplicacion(int dimension, double a[dimension][dimension], double b[dimension][dimension], double c[dimension][dimension]);

void transpuesta(int dimension, double a[dimension][dimension],  double transpuesta[dimension][dimension]);
#endif // MATLABUTILITIES_H_INCLUDED
