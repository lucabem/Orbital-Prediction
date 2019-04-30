
#include "Constantes.h"
#include <math.h>


void Position (long double lon, long double lat, long double h, long double pos[3])
{
    long double R_equ = R_Earth;
    long double f     = f_Earth;



    long double e2     = f*(2-f);
    // Square of eccentricity
    long double CosLat = cos(lat);
    // (Co)sine of geodetic latitude
    long double SinLat = sin(lat);

    // Position vector
    long double N = R_equ/sqrt(1-e2*SinLat*SinLat);

    pos[0] =  (N+h)*CosLat*cos(lon);
    pos[1] =  (N+h)*CosLat*sin(lon);
    pos[2] =  ((1-e2)*N+h)*SinLat;


}
