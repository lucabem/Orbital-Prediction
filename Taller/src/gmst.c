
#include "..\gmst.h"
#include "..\Frac.h"
#include <math.h>
#include "..\Constantes.h"

/**
    Funcion que calcula el tiempo Greenwich sideral medio.

    Entrada:
        double Mjd_UT1:  Fecha juliana modificada

    Salida:
        double gmstime:  GMST en [rad]

*/
 double gmst( double Mjd_UT1)
{
    int secs = 86400;

     double Mjd_0 = floor(Mjd_UT1);
     double UT1   = secs*(Mjd_UT1-Mjd_0);

     double T_0   = (Mjd_0  - MJD_J2000)/36525;
     double T     = (Mjd_UT1 - MJD_J2000)/36525;

     double gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1  + (0.093104-6.2e-6*T)*T*T;

     double gmstime = 2*PI*Frac(gmst/secs);

    return gmstime;

}
