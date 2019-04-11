#include "MatlabUtilities.h"


double Norma(double v[])
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

double dot(double v1[], double v2[])
{
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
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

void multiplicacion(double a[3][3], double b[3][3], double c[3][3])
{
    int i,j,k;
    for(i=0; i<=3-1; i++)
    {
        for(j=0; j<=3-1; j++)
        {
            for(k=0; k<=3-1; k++)
            {
                c[i][j]=(c[i][j]+(a[i][k]*b[k][j]));
            }
        }
    }
    return c;
}


