
#include <math.h>
#include "Frac.h"

/**
    Funcion que dado un x devuelve su parte fracionaria.

    Entrada:
        x - double
    Salida
        Parte fraccionaria

*/
 double Frac(  double x)
{
    return x-floor(x);
}
