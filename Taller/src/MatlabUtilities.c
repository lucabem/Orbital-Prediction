#include "..\MatlabUtilities.h"


/**
    Funcion que calcular la norma de un vector de dimension 3.
*/
 double Norma( double v[])
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}


/**
    Funcion que calcula el producto escalar entre dos vectores de la misma dimension.
*/
 double dot(int dim,  double v1[dim],  double v2[dim])
{
     double suma = 0;
    for(int i=0; i<dim; i++)
    {
        suma = suma + v1[i]*v2[i];
    }
    return suma;
}

/**
    Funcion que devuelve el signo de x:
         1 si x > 0
        -1 si x < 0
         0 si x = 0
*/
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


/**
    Accion que calcula el vector cr3, resultado del producto vectorial entre dos vectores de la misma dimension.
*/
void cross( double v1[],  double v2[],  double cr[3])
{
    cr[0] = v1[1]*v2[2]-v1[2]*v2[1];
    cr[1] = v1[2]*v2[0]-v1[0]*v2[2];
    cr[2] = v1[0]*v2[1]-v1[1]*v2[0];

}

/**

    Devuelve el determinante de una matriz cuadrada. Solo se permiten matrices de dimension 3x3 y 2x2.
    En caso de otra dimension, devuelve 0.

*/

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


/**
    Acción que multiplica la matriz A de dimension r1 x c1 y la matriz B de dimension r2 x c2.
    Guarda el resultado en la matriz C de dimension r1 x c2.
*/
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

/**
    Accion que dado el numero de filas, columnas y la matriz ha transponer, genera la matriz transpuesta
    y la guarda en el parametro transpuesta.
*/
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

/**
    Funcion que comprueba si un número es real o imaginario.
    @return True si parteImaginaria es distinto de 0. False si parteImaginaria es igual a 0.
*/
bool isReal( double parteImaginaria)
{
    return (parteImaginaria) == 0;
}
/**
    Accion que dada una matrz de dimensión M x N,
    inicializa todos sus componentes a cero.
*/

void zeros(int m, int n,  double matriz[m][n])
{
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            matriz[i][j] = 0.0;
}


/**

    Funcion que calcular las raices (complejas y reales) de un polinomio de grado n y las guarda en el vector resultado.

    @param[degree]: grado del polinomio
    @param[op]: vector con los coeficientes del polinomio
    @param[zeror]: vector con las partes reales de las soluciones
    @param[zeroi]: vector con las partes imaginarias de las soluciones
    @param[resultado]: vector con las raices reales del polinomio


*/
int raicesPolinomiales(int degree,  double op[degree+1], double zeror[degree+1],  double zeroi[degree+1], double resultado[20])
{

    int info[15] = {1,2,3,4,5,6,7, 8, 9, 10, 11, 12, 13, 14, 15};

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


/**
    Funcion que dado un double x:
        Si x > 0, devuelve la parte entera de x.
        Si x < 0, devuelve la parte entera de x + 1
*/
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


/**
    Funcion para comparar dos reales con una precision de 10^(-12).
    @return True si  abs(a-b) es menor que 10^(-12), y False en caso contrario.
*/
bool fequal( double a,  double b)
{
    return fabs(a-b) < 0.0000000000001;
}

double modulo(double a, double b){ return (((a/b) - floor(a/b))*b ); }

