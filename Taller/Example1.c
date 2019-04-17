
#include "Constantes.h"
#include <stdio.h>
#include <stdlib.h>
#include "MatlabUtilities.h"

double Mjday(int year, int month,int day, int hour, int min, double sec);
void Position (double lon, double lat, double h, double pos[3]);


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
            i++;
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
        obs[i][1] = Rad*rtasc;
        obs[i][2] = Rad*decl;
        i = i+1;

    }

    fclose(fid);

    double lat = Rad*39.13607;     // [rad]
    double lon = Rad*(-121.35072); // [rad]
    double alt = 0.09981638e3;     // [m]

    double Rs[3];
    Position(lon, lat, alt, Rs);

    double Mjd1 = obs[0][0];
    double Mjd2 = obs[1][0];
    double Mjd3 = obs[2][0];

    double Mjd_UTC = Mjd1;

    return 0;

}


double Mjday(int year, int month,int day, int hour, int min, double sec)
{

    int y = year;
    int m = month;
    double b = 0.0;
    double c = 0.0;
    double a = 0.0;
    if (m <=2)
    {
        y = y-1;
        m = m+12;
    }
    if(y<0)
    {
        c = -0.75;
    }

    // check for valid calendar date
    if (year < 1582)
    {

    }
    else if (year > 1582)
    {
        a = fix(y/100);
        b = 2 - a + floor(a/4);

    }
    else if ( month < 10)
    {

    }
    else if (month > 10)
    {
        a = fix(y/100);
        b = 2 - a + floor(a/4);
    }
    else if (day <= 4)
    {

    }
    else if (day > 14)
    {
        a = fix(y/100);
        b = 2 - a + floor(a/4);
    }
    else
    {
        return 0.0;
    }


    double jd = fix(365.25*y+c) + fix(30.6001 * (m + 1));

    jd = jd + day + b + 1720994.5;

    jd = jd + (hour*1.0+min*1.0/60+sec/3600)/24;

    double Mjd = jd - 2400000.5;

    return Mjd;
}


void Position (double lon, double lat, double h, double pos[3])
{
    double R_equ = R_Earth;
    double f     = f_Earth;

    double e2     = f*(2-f);
    // Square of eccentricity
    double CosLat = cos(lat);
    // (Co)sine of geodetic latitude
    double SinLat = sin(lat);

    // Position vector
    double N = R_equ/sqrt(1-e2*SinLat*SinLat);

    pos[0] =  (N+h)*CosLat*cos(lon);
    pos[1] =  (N+h)*CosLat*sin(lon);
    pos[2] =  ((1-e2)*N+h)*SinLat;

}

