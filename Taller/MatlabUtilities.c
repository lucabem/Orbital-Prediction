#include "MatlabUtilities.h"

double Norma(double v[])
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

double dot(int dim, double v1[dim], double v2[dim])
{
    double suma = 0;
    for(int i=0; i<dim; i++)
    {
        suma = suma + v1[i]*v2[i];
    }
    return suma;
}

int Sing(double x)
{
    if ( x<0 )
    {
        return -1;
    }
    else if( x>0 )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

double* cross(double v1[], double v2[])
{

    double static resultado[3] = {0.0, 0.0, 0.0};
    resultado[0] = v1[1]*v2[2]-v1[2]*v2[1];
    resultado[1] = v1[2]*v2[0]-v1[0]*v2[2];
    resultado[2] = v1[0]*v2[1]-v1[1]*v2[0];

    return resultado;
}

double det(int dimension, double matrix[dimension][dimension])
{
    double det = 0;
    if(dimension == 3)
    {
        for(int i=0; i<3; i++)
            det = det + (matrix[0][i]*(matrix[1][(i+1)%3]*matrix[2][(i+2)%3] - matrix[1][(i+2)%3]*matrix[2][(i+1)%3]));
    }

    if(dimension == 2)
    {
        det = matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0];
    }

    return det;
}

void multiplicacion(int r1, int c1, int r2, int c2, double a[r1][c1], double b[r2][c2], double c[r1][c2])
{
    int i,j,k;
    for(i=0; i<r1; ++i)
    {
        for(j=0; j<c2; ++j)
        {
            for(k=0; k<c1; ++k)
            {
                c[i][j]+=a[i][k]*b[k][j];
            }
        }
    }

}

void transpuesta(int filas, int columnas, double a[filas][columnas], double transpuesta[columnas][filas])
{

    for(int i=0; i<filas; i++)
    {
        for(int j=0; j<columnas; j++)
        {
            transpuesta[j][i] = a[i][j];
        }
    }

}

bool isReal(double parteImaginaria)
{
    return (parteImaginaria) == 0;
}

double* zeros(int m, int n)
{
    double *mat;
    mat = calloc(m*n, sizeof(double));
    return mat;
}

double* roots(double *op, int degree, double *zeror, double *zeroi)
{

    int info[7] = {1,2,3,4,5,6,7};
    int contado = 0;

    rpoly(op,degree, zeror, zeroi, info);

    for(int i=0; i<degree; i++)
    {
        if (isReal(zeroi[i]))
        {
            contado++;
        }
    }

    double resultado[contado];
    int posicion = 0;
    for(int i=0; i<degree; i++)
    {
        if (isReal(zeroi[i]))
        {
            resultado[posicion] = zeror[i];
            posicion++;
        }
    }

    return resultado;

}

double fix(double x)
{
    if (x >= 0)
    {
        return floor(x);
    }
    else
    {
        return ceil(x);
    }

}
