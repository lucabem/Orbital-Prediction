#include "MatlabUtilities.h"

 double Norma( double v[])
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

 double dot(int dim,  double v1[dim],  double v2[dim])
{
     double suma = 0;
    for(int i=0; i<dim; i++)
    {
        suma = suma + v1[i]*v2[i];
    }
    return suma;
}

int Sing( double x)
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

void cross( double v1[],  double v2[],  double cr[3])
{
    cr[0] = v1[1]*v2[2]-v1[2]*v2[1];
    cr[1] = v1[2]*v2[0]-v1[0]*v2[2];
    cr[2] = v1[0]*v2[1]-v1[1]*v2[0];

}

 double det(int dimension,  double matrix[dimension][dimension])
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

void multiplicacion(int r1, int c1, int r2, int c2,  double a[r1][c1],  double b[r2][c2],  double c[r1][c2])
{
    int i,j,k;
    for(int x=0; x<r1; x++)
        for(int y=0; y<c2; y++)
            c[x][y] = 0.0;


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

void transpuesta(int filas, int columnas,  double a[filas][columnas],  double transpuesta[columnas][filas])
{
    for(int i=0; i<filas; i++)
    {
        for(int j=0; j<columnas; j++)
        {
            transpuesta[j][i] = a[i][j];
        }
    }
}

bool isReal( double parteImaginaria)
{
    return (parteImaginaria) == 0;
}

void zeros(int m, int n,  double matriz[m][n])
{
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            matriz[i][j] = 0.0;
}

int raicesPolinomiales(int degree,  double op[degree+1], double zeror[degree+1],  double zeroi[degree+1], double resultado[20])
{

    int info[15] = {1,2,3,4,5,6,7, 8, 9, 10, 11, 12, 13, 14, 15};
    int contado = 0;

    rpoly(op,degree, zeror, zeroi, info);

    int posicion = 0;

    for(int i=0; i<degree; i++)
    {
        if (isReal(zeroi[i]))
        {
            resultado[posicion] = zeror[i];
            posicion++;
        }
    }

    return 0;
}

 double fix( double x)
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

bool fequal( double a,  double b)
{
    return fabs(a-b) < 0.0000000000001;
}

