
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
        eop[fila][0] = v1; eop[fila][1] = v2; eop[fila][2] = v3; eop[fila][3] = v4; eop[fila][4] = v5; eop[fila][5] = v6;
        eop[fila][6] = v7; eop[fila][7] = v8; eop[fila][8] = v9; eop[fila][9] = v10; eop[fila][10] = v11; eop[fila][11] = v12;
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

    double P[3][3], N[3][3], PoleM[3][3];

    PrecMatrix(MJD_J2000, Mjd_TT, P);

    NutMatrix(Mjd_TT, N);

    PoleMatrix(x_pole, y_pole, PoleM);

    free(eop);
    return 0;

}






