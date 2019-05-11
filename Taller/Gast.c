#include "..\Gast.h"
#include "..\timeDiff.h"
#include "..\gmst.h"
#include "..\EqnEquinox.h"
#include <math.h>
#include "..\IERS.h"
#include "..\Constantes.h"

/**
    Funcion que a partir de la fecha juliana modificada calcula
    el tiempo sideral Greenwich aparente.

    Entrada:
        Mjd_UT1 Fecha Julia Modificada

    Salida:
        GAST en radianes

*/
 long double gast ( long double Mjd_UT1,  long double(*eop)[13])
{
    long double salidaIERS[6], diferenciaTiempos[5];

    IERS(eop, Mjd_UT1, 'l', salidaIERS);

    timeDiff(salidaIERS[0], salidaIERS[1], diferenciaTiempos);

     long double Mjd_UTC = Mjd_UT1 - salidaIERS[0]/86400;
     long double Mjd_TT = Mjd_UTC + diferenciaTiempos[3]/86400;


    return fmod(gmst(Mjd_UT1) + eqnEquinox(Mjd_TT), PI2);
}
