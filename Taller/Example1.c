
#include "Constantes.h"
#include <stdio.h>
#include <stdlib.h>
#include "MatlabUtilities.h"
#include "Mjday.h"
#include "PrecMatrix.h"
#include "Position.h"

int main()
{
    FILE* fid = fopen("eop19620101.txt","rt");

    int v1, v2, v3, v4, v13;
    float v5, v6, v7, v8, v9, v10, v11, v12;
    if (fid == NULL)
    {
        exit(EXIT_FAILURE);
    }
    /*
        while( fscanf(fid,"%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12, &v13) != EOF)
        {
                printf("%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d \n", v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13);
        }
    */
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


    return 0;

}





