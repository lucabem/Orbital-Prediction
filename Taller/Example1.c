
#include "Constantes.h"
#include <stdio.h>
#include <stdlib.h>
#include "MatlabUtilities.h"
#include "Mjday.h"
#include "PrecMatrix.h"
#include "Position.h"
#include "rv2coe.h"
#include "NutAngles.h"
#include "NutMatrix.h"
#include "EqnEquinox.h"
#include "Gibbs.h"
#include "NewtonNu.h"

/*

int main()
{
    FILE* fid = fopen("eop19620101.txt","rt");

    int v1, v2, v3, v4, v13;
    float v5, v6, v7, v8, v9, v10, v11, v12;
    if (fid == NULL)
    {
        exit(EXIT_FAILURE);
    }

        while( fscanf(fid,"%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12, &v13) != EOF)
        {
                printf("%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d \n", v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13);
        }

    fclose(fid);

    fid = fopen("sat1.txt", "rt");
    int i = 0;
    int Y, M, D, h, m;
    float s, rtasc, decl;

    double obs[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    while(1)
    {
        if (feof(fid))
        {
            break;
        }

        fscanf(fid, "%d/%d/%d %d:%d:%f  %f   %f", &Y, &M, &D, &h, &m, &s, &rtasc, &decl);

        obs[i][0] = Mjday(Y,M,D,h,m,s);
        obs[i][1] = RAD*rtasc;
        obs[i][2] = RAD*decl;
        i = i+1;

    }

    fclose(fid);

    double lat = RAD*39.13607;     // [rad]
    double lon = RAD*(-121.35072); // [rad]
    double alt = 0.09981638e3;     // [m]

    double Rs[3];
    Position(lon, lat, alt, Rs); // vector [a b c] -> matrix 1x3

    double Mjd1 = obs[0][0];
    double Mjd2 = obs[1][0];
    double Mjd3 = obs[2][0];

    double Mjd_UTC = Mjd1;

    double matrizTran[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};


    PrecMatrix(MJD_J2000, 54977.66766966425, matrizTran);


    double r[3] = {20435422.3521544, 1070699.44671825, 1012905.49143365};
    double v[3] = {17.1964697862374, -2657.51027611478, 3738.38685080781};
    double resultado[12] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    rv2coe(r, v, resultado);

    double r[2] = {0.0, 0.0};
    NutAngles(54977.667670, r);
    printf("Salida  dpsi=%f  deps=%f \n", r[0], r[1]);


    NutMatrix(54977.6746141088, matrizTran);

    for(int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            printf("%0.15f ", matrizTran[i][j]);
        }
        printf("\n");
    }



    printf("Resultado -> %0.17f  \n", eqnEquinox(54977.6676696643));

    double r1[3] = {  20387627.0717529,
                      1865163.69633398,
                      -109943.688555879
                   };
    double r2[3] = { 20435422.3521544,
                     1070699.44671825,
                     1012905.49143365
                   };
    double r3[3] = { 20398157.0666256,
                     271778.615869788,
                     2131538.39542076
                   };

    double salida[3];

    char *error = "";
    double angTheta[2];
    gibbs(r1, r2, r3, salida, angTheta, error);

    printf("theta = %0.15f \n", angTheta[0]);

    for (int i=0; i<3; i++)
        printf("%0.15f \n", salida[i]);



    double salida[2];
    newtonnu(0.0825331061733743, 0.181200311069884, salida);

    for(int i=0; i<2; i++)
        printf("%0.15f \n", salida[i]);


    return 0;

}

*/




