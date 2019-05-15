
#include "Constantes.h"
#include <math.h>

/**

    Posici�n: vector de posici�n (r [m]) de coordenadas geod�sicas
             (Longitud [rad], latitud [rad], altitud [m])

    Last modified:   2015/08/12   M. Mahooti
*/

void Position ( double lon,  double lat,  double h,  double pos[3])
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
