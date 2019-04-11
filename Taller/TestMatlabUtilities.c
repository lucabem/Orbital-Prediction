#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "MatlabUtilities.h"
#include <math.h>

void Test_Norma();
void Test_Sign();
void Test_Cross();
void Test_Dot();
void Test_Det();
void Test_Multiplicacion();
void Test_Transpuesta();

bool fequal(double a, double b);

const double epsilon = 0.00001;

int main()
{
    Test_Norma();
    Test_Sign();
    Test_Cross();
    Test_Dot();
    Test_Det();
    Test_Multiplicacion();
    Test_Transpuesta();

    printf(">>> TESTS MATLABUTILITIES SUPERADOS \n");
    return 0;
}

void Test_Norma()
{
    double v1[3] = {0.0, 0.0, 0.0};
    double r1 = Norma(v1);
    assert(fequal(0.0, r1) == true);

    double v2[3] = {1, 1, 1};
    double r2 = Norma(v2);
    assert(fequal(1.73205, r2) == true);

    double v3[3] = {4.0, 3.0, 0.0};
    double r3 = Norma(v3);
    assert(fequal(5.0, r3) == true);
}

void Test_Sign()
{
    assert(Sing(0.0) == 0);
    assert(Sing(1.5) == 1);
    assert(Sing(-0.5) == -1);
}

void Test_Cross()
{

    double v1[3] = {4.0, -2.0, 1.0};
    double v2[3] = {1.0, -1.0, 3.0};

    double* resultado = cross(v1, v2);

    assert(fequal(resultado[0], -5.00) == true);
    assert(fequal(resultado[1], -11.00) == true);
    assert(fequal(resultado[2], -2.00) == true);

    double v3[3] = {0.25, -2.75, 12.2223};
    double v4[3] = {1.045, -1.12, 3.5};

    resultado = cross(v3, v4);

    assert(fequal(resultado[0], 4.063976) == true);
    assert(fequal(resultado[1], 11.8973035) == true);
    assert(fequal(resultado[2], 2.59375) == true);
}

void Test_Dot()
{
    double v1[3] = {4.0, -1.0, 2.0};
    double v2[3] = {2.0, -2.0, -1.0};
    assert(fequal(dot(v1, v2), 8) == true);

    double v3[3] = {10.1220, 0.0, -12.34560};
    double v4[3] = {-5.35, 2.67, 0.5};
    assert(fequal(dot(v3, v4), -60.3255) == true);
}

void Test_Det()
{
    double matrix[3][3] = {{1, -2, 4}, {-5,2,0}, {1,0,3}};
    double determinante = det(3, matrix);

    assert(fequal(determinante, -32 ) == true);

    double matrix2[3][3] = {{1.056, -2.778, 4.123}, {-5.35, 2.67, 0.5}, {1.11, 0, 3.23}};
    double determinante2 = det(3, matrix2);

    assert(fequal(determinante2, -52.6593 ) == true);

    double matrix3[2][2] = {{1.056, -2.778}, {-5.35, 2.67}};
    double determinante3 = det(2, matrix3);

    assert(fequal(determinante3, -12.04278 ) == true);


}

void Test_Multiplicacion()
{
    double a[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    double c[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    multiplicacion(3, a, a, c);

    for(int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            if (i==j)
            {
                assert(fequal(c[i][j], 1.0) == true);
            }
            else
            {
                assert(fequal(c[i][j], 0.0) == true);
            }
        }
    }

    double b[3][3] = {{1.056, -2.778, 4.123}, {-5.35, 2.67, 0.5}, {1.11, 0, 3.23}};
    double d[3][3] = {{10.1220, 0.0, -12.34560}, {0.25, -2.75, 12.2223}, {-5.35, 2.67, 0.5}};
    double f[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    multiplicacion(3, b,d, f);

    assert(fequal(f[0][0], -12.063718) == true);
    assert(fequal(f[0][1], 18.64791) == true);
    assert(fequal(f[0][2], -44.929003) == true);
    assert(fequal(f[1][0], -56.1602) == true);
    assert(fequal(f[1][1], -6.0075) == true);
    assert(fequal(f[1][2], 98.932501) == true);
    assert(fequal(f[2][0], -6.04508) == true);
    assert(fequal(f[2][1], 8.6241) == true);
    assert(fequal(f[2][2], -12.088616) == true);

    double x[2][2] = {{1.056, -2.778}, {-5.35, 2.67}};
    double y[2][2] = {{10.1220, -12.34560}, {0.25, -2.75}};
    double z[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

    multiplicacion(2, x, y, z);
    assert(fequal(z[0][0], 9.994332) == true);
    assert(fequal(z[0][1], -5.3974536) == true);
    assert(fequal(z[1][0], -53.4852) == true);
    assert(fequal(z[1][1], 58.70646) == true);

}

void Test_Transpuesta()
{

    double matrix[3][3] = {{1, -2, 4}, {-5,2,0}, {1,0,3}};
    double trans[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    transpuesta(3,matrix, trans);

    assert(fequal(trans[0][0], 1) == true);
    assert(fequal(trans[0][1], -5) == true);
    assert(fequal(trans[0][2], 1) == true);
    assert(fequal(trans[1][0], -2) == true);
    assert(fequal(trans[1][1], 2) == true);
    assert(fequal(trans[1][2], 0) == true);
    assert(fequal(trans[2][0], 4) == true);
    assert(fequal(trans[2][1], 0) == true);
    assert(fequal(trans[2][2], 3) == true);

    double mat[2][2] = {{1.24, 2.7800}, {0.0050, 1.1120}};
    double t[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

    transpuesta(2, mat, t);

    assert(fequal(t[0][0], 1.24) == true);
    assert(fequal(t[0][1], 0.005) == true);
    assert(fequal(t[1][0], 2.78) == true);
    assert(fequal(t[1][1], 1.112) == true);

}

bool fequal(double a, double b)
{
    return fabs(a-b) < epsilon;
}



