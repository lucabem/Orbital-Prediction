
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


int main()
{

    double (*eop)[13] = malloc(sizeof(double[20026][13]));

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

//    assert( fequal(0.999997373802329, P[0][0]) == true);
//    assert( fequal(-0.00210194819368335, P[0][1]) == true);
//    assert( fequal(-0.00091334672251536, P[0][2]) == true);
//    assert( fequal(0.00210194819366918, P[1][0]) == true);
//    assert( fequal(0.999997790903995, P[1][1]) == true);
//    assert( fequal(-0.000000959920515567498, P[1][2]) == true);
//    assert( fequal(0.000913346722547957, P[2][0]) == true);
//    assert( fequal(-0.000000959889498958352, P[2][1]) == true);
//    assert( fequal(0.999999582898335, P[2][2]) == true);
//
//
//    assert( fequal(0.999999997895984, N[0][0]) == true);
//    assert( fequal(-0.0000595170051004573, N[0][1]) == true);
//    assert( fequal(-0.0000258022713429399, N[0][2]) == true);
//    assert( fequal(0.0000595164295625406, N[1][0]) == true);
//    assert( fequal(0.999999997980121, N[1][1]) == true);
//    assert( fequal(-0.0000223059014984558, N[1][2]) == true);
//    assert( fequal(0.0000258035988712757, N[2][0]) == true);
//    assert( fequal(0.0000223043657924436, N[2][1]) == true);
//    assert( fequal(0.999999999418345, N[2][2]) == true);
//    printf("\n");
//
//
//    assert( fequal(1.0, PoleM[0][0]) == true);
//    assert( fequal(0.000000000000194609382460874, PoleM[0][1]) == true);
//    assert( fequal(0.0000000757892008065431, PoleM[0][2]) == true);
//    assert( fequal(0.0, PoleM[1][0]) == true);
//    assert( fequal(0.999999999996703, PoleM[1][1]) == true);
//    assert( fequal(-0.00000256777193042298, PoleM[1][2]) == true);
//    assert( fequal(-0.0000000757892008067929, PoleM[2][0]) == true);
//    assert( fequal(0.00000256777193042298, PoleM[2][1]) == true);
//    assert( fequal(0.9999999999967, PoleM[2][2]) == true);
//    printf("\n");

//  PoleM*GHA*N*P;
//  --------------
//  resultadoPG = PoleM*GHA
//  resultadoPGN = resultadoPG*N
//  E = resultadoPGN*P

    multiplicacion(3, 3, 3, 3, PoleM, GHA, resultadoPG);
    multiplicacion(3, 3, 3, 3, resultadoPG, N, resultadoPGN);
    multiplicacion(3, 3, 3, 3, resultadoPGN, P, E);



//     assert( fequal(-0.563792466027701, E[0][0]) == true);
//     assert( fequal(0.825916335412934, E[0][1]) == true);
//     assert( fequal(0.000512004369342484, E[0][2]) == true);
//     assert( fequal(-0.825915962467244, E[1][0]) == true);
//     assert( fequal(-0.563792698167355, E[1][1]) == true);
//     assert( fequal(0.000785133734051811, E[1][2]) == true);
//     assert( fequal(0.000937119101302219, E[2][0]) == true);
//     assert( fequal(0.0000197799025896599, E[2][1]) == true);
//     assert( fequal(0.999999560708176, E[2][2]) == true);

    printf("\n");
    transpuesta(3, 3, E, Et);

    double RsMatriz[3][1] = {{Rs[0]}, {Rs[1]}, {Rs[2]}};
    double rsite1Matriz[3][1];
    multiplicacion(3, 3, 3, 1, Et, RsMatriz, rsite1Matriz);

    double rsite1[3] = {rsite1Matriz[0][0], rsite1Matriz[1][0], rsite1Matriz[2][0]};

    printf("\n ----- CALCULO RSITE1 ----- \n");
    for (int i=0; i<3; i++)
        printf("%0.12f \n", rsite1[i]);
    printf("\n -------------------- \n");

/*
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

        printf("\n -----CALCULO RSITE2----- \n");
        for (int i=0; i<3; i++)
            printf("%0.5f \n", rsite2[i]);
        printf("\n -------------------- \n");

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

        printf("\n -----CALCULO RSITE3----- \n");
        for (int i=0; i<3; i++)
            printf("%0.5f \n", rsite3[i]);
        printf("\n -------------------- \n");

*/

    free(eop);
    return 0;

}






