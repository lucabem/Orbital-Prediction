
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
#include "HGibbs.h"
#include "DoubleR.h"
#include "IERS.h"
#include "timeDiff.h"
#include "PoleMatrix.h"
#include "Gast.h"
#include "GHAMatrix.h"
#include <assert.h>
#include "lambert_gooding.h"
#include "AngleDr.h"


void Example6();


void Example6()
{

    double (*eop)[13] = malloc(sizeof( double[20026][13]));

    FILE* fid = fopen("eop19620101.txt","rt");

    int v1, v2, v3, v4, v13;
    float v5, v6, v7, v8, v9, v10, v11, v12;
    if (fid == NULL)
    {
        exit(EXIT_FAILURE);
    }

    int fila = 0;
    while( fscanf(fid,"%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12, &v13) != EOF)
    {
        eop[fila][0] = v1;
        eop[fila][1] = v2;
        eop[fila][2] = v3;
        eop[fila][3] = v4;
        eop[fila][4] = v5;
        eop[fila][5] = v6;
        eop[fila][6] = v7;
        eop[fila][7] = v8;
        eop[fila][8] = v9;
        eop[fila][9] = v10;
        eop[fila][10] = v11;
        eop[fila][11] = v12;
        eop[fila][12] = v13;

        fila++;
    }
    fclose(fid);

    fid = fopen("sat6.txt", "rt");
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

        fscanf(fid, "%d/%d/%d %d:%d:%f %f %f", &Y, &M, &D, &h, &m, &s, &rtasc, &decl);

        obs[i][0] = Mjday(Y,M,D,h,m,s);
        obs[i][1] = RAD*rtasc;
        obs[i][2] = RAD*decl;
        i = i+1;

    }

    fclose(fid);

    double lat = RAD*30.5724;     // [rad]
    double lon = RAD*(-86.2143); // [rad]
    double alt = 0.0;     // [m]


    double Rs[3];
    Position(lon, lat, alt, Rs); // vector [a b c] -> matrix 1x3

    double Mjd1 = obs[0][0];
    double Mjd2 = obs[1][0];
    double Mjd3 = obs[2][0];

    double Mjd_UTC = Mjd1;

    double salida[6];
    double UT1_UTC, TAI_UTC, x_pole, y_pole, ddpsi, ddeps;
    IERS(eop, Mjd_UTC, 'l', salida);
    UT1_UTC = salida[0];
    TAI_UTC = salida[1];
    x_pole  = salida[2];
    y_pole  = salida[3];
    ddpsi   = salida[4];
    ddeps   = salida[5];

    double diferenciaTiempos[5];

    timeDiff(UT1_UTC, TAI_UTC, diferenciaTiempos);

    double UT1_TAI = diferenciaTiempos[0];
    double UTC_GPS = diferenciaTiempos[1];
    double UT1_GPS = diferenciaTiempos[2];
    double TT_UTC  = diferenciaTiempos[3];
    double GPS_UTC = diferenciaTiempos[4];

    double Mjd_TT = Mjd_UTC + TT_UTC/86400;
    double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    double P[3][3], N[3][3], PoleM[3][3], GHA[3][3];

    PrecMatrix(MJD_J2000, Mjd_TT, P);

    NutMatrix(Mjd_TT, N);

    PoleMatrix(x_pole, y_pole, PoleM);

    ghaMatrix(Mjd_UT1, GHA, eop);

    double resultadoPG[3][3], resultadoPGN[3][3], E[3][3], Et[3][3];

    multiplicacion(3, 3, 3, 3, PoleM, GHA, resultadoPG);
    multiplicacion(3, 3, 3, 3, resultadoPG, N, resultadoPGN);
    multiplicacion(3, 3, 3, 3, resultadoPGN, P, E);

    transpuesta(3, 3, E, Et);

    double RsMatriz[3][1] = {{Rs[0]}, {Rs[1]}, {Rs[2]}};
    double rsite1Matriz[3][1];
    multiplicacion(3, 3, 3, 1, Et, RsMatriz, rsite1Matriz);

    double rsite1[3] = {rsite1Matriz[0][0], rsite1Matriz[1][0], rsite1Matriz[2][0]};

    Mjd_UTC = Mjd2;

    IERS(eop, Mjd_UTC, 'l', salida);
    UT1_UTC = salida[0];
    TAI_UTC = salida[1];
    x_pole  = salida[2];
    y_pole  = salida[3];
    ddpsi   = salida[4];
    ddeps   = salida[5];


    timeDiff(UT1_UTC, TAI_UTC, diferenciaTiempos);

    UT1_TAI = diferenciaTiempos[0];
    UTC_GPS = diferenciaTiempos[1];
    UT1_GPS = diferenciaTiempos[2];
    TT_UTC  = diferenciaTiempos[3];
    GPS_UTC = diferenciaTiempos[4];

    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;


    PrecMatrix(MJD_J2000, Mjd_TT, P);

    NutMatrix(Mjd_TT, N);

    PoleMatrix(x_pole, y_pole, PoleM);

    ghaMatrix(Mjd_UT1, GHA, eop);

    multiplicacion(3, 3, 3, 3, PoleM, GHA, resultadoPG);
    multiplicacion(3, 3, 3, 3, resultadoPG, N, resultadoPGN);
    multiplicacion(3, 3, 3, 3, resultadoPGN, P, E);

    transpuesta(3, 3, E, Et);

    double RsMatriz2[3][1] = {{Rs[0]}, {Rs[1]}, {Rs[2]}};
    double rsite2Matriz[3][1];
    multiplicacion(3, 3, 3, 1, Et, RsMatriz2, rsite2Matriz);

    double rsite2[3] = {rsite2Matriz[0][0], rsite2Matriz[1][0], rsite2Matriz[2][0]};


    Mjd_UTC = Mjd3;

    IERS(eop, Mjd_UTC, 'l', salida);
    UT1_UTC = salida[0];
    TAI_UTC = salida[1];
    x_pole  = salida[2];
    y_pole  = salida[3];
    ddpsi   = salida[4];
    ddeps   = salida[5];


    timeDiff(UT1_UTC, TAI_UTC, diferenciaTiempos);

    UT1_TAI = diferenciaTiempos[0];
    UTC_GPS = diferenciaTiempos[1];
    UT1_GPS = diferenciaTiempos[2];
    TT_UTC  = diferenciaTiempos[3];
    GPS_UTC = diferenciaTiempos[4];

    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;


    PrecMatrix(MJD_J2000, Mjd_TT, P);

    NutMatrix(Mjd_TT, N);

    PoleMatrix(x_pole, y_pole, PoleM);

    ghaMatrix(Mjd_UT1, GHA, eop);

    multiplicacion(3, 3, 3, 3, PoleM, GHA, resultadoPG);
    multiplicacion(3, 3, 3, 3, resultadoPG, N, resultadoPGN);
    multiplicacion(3, 3, 3, 3, resultadoPGN, P, E);

    transpuesta(3, 3, E, Et);

    double RsMatriz3[3][1] = {{Rs[0]}, {Rs[1]}, {Rs[2]}};
    double rsite3Matriz[3][1];
    multiplicacion(3, 3, 3, 1, Et, RsMatriz3, rsite3Matriz);

    double rsite3[3] = {rsite3Matriz[0][0], rsite3Matriz[1][0], rsite3Matriz[2][0]};

    double vv2[3], r2[3];

    anglesdr(obs[0][1], obs[1][1], obs[2][1], obs[0][2], obs[1][2], obs[2][2],
             Mjd1, Mjd2, Mjd3, rsite1, rsite2, rsite3, r2, vv2);

    double errorVV2[3] = { 0.607488210306698, -1.75441167132057 , 6.58484848503918};
    double errorR2[3] = {  8261.9716892784, 2253.75362052551 , 1969.11103974355};

/*
    printf("\n RESULTADO FINAL EXAMPLE 6 Double-R-Iteration method\n");
    printf("\t%0.5f ------- %0.5f \n ", r2[0]/1000, errorR2[0]);
    printf("\t%0.5f  ------- %0.5f \n ", r2[1]/1000, errorR2[1]);
    printf("\t%0.5f ------- %0.5f \n ", r2[2]/1000, errorR2[2]);
    printf("\t%0.5f    ------- %0.5f \n ", vv2[0]/1000, errorVV2[0]);
    printf("\t%0.5f    ------- %0.5f \n ", vv2[1]/1000, errorVV2[1]);
    printf("\t%0.5f   ------- %0.5f \n ", vv2[2]/1000, errorVV2[2]);
*/
    printf("ERROR APROXIMADO \n");
    printf("\t%0.5f \n ", fabs(r2[0]/1000 - errorR2[0])*1000);
    printf("\t%0.5f \n ", fabs(r2[1]/1000 - errorR2[1])*1000);
    printf("\t%0.5f \n ", fabs(r2[2]/1000 - errorR2[2])*1000);
    printf("\t%0.5f \n ", fabs(vv2[0]/1000 - errorVV2[0])*1000);
    printf("\t%0.5f \n ", fabs(vv2[1]/1000 - errorVV2[1])*1000);
    printf("\t%0.5f \n ", fabs(vv2[2]/1000 - errorVV2[2])*1000);

    free(eop);
}
