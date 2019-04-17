
#include "Constantes.h"
#include <math.h>
#include "MatlabUtilities.h"

void IERS(double* eop, double Mjd_UTC, char interp)
{
    double Arcs = 3600*18/PI;

    if (interp == '1')
    {
        double mj = floor(Mjd_UTC);
        int nop = NELEMS(eop)/13; //20026

        double preeop, nexteop;

    }

}
