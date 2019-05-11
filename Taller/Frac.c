
#include <math.h>
#include "Frac.h"

/**
    Funcion que dado un x devuelve su parte fracionaria.

    Entrada:
        x - long double
    Salida
        Parte fraccionaria

*/
 long double Frac(  long double x)
{
    return x-floor(x);
}
