
#include "Constantes.h"
#include <math.h>
#include "MatlabUtilities.h"
#include "IERS.h"

void IERS(double(*eop)[13], double Mjd_UTC, char interp, double salida[6])
{
    double Arcs = 3600*18/PI;

    if (interp == 'l')
    {
        double mj = floor(Mjd_UTC);
        int nop = 20026;

        double preeop[13], nexteop[13];

        for (int i=0; i<nop; i++)
        {
            if (mj == eop[i][3])
            {
                for (int j=0; j<13; j++)
                {
                    preeop[j] = eop[i][j] * 1.0;
                    nexteop[j] = eop[i+1][j] * 1.0;
                }
                break;
            }
        }

        double mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
        double fixf = mfme/1440;

        double UT1_UTC = preeop[6] + (nexteop[6]-preeop[6])*fixf;
        double TAI_UTC = preeop[12];
        double x_pole  = preeop[4]+(nexteop[4]-preeop[4])*fixf;
        double y_pole  = preeop[5]+(nexteop[5]-preeop[5])*fixf;
        double ddpsi   = preeop[8]+(nexteop[8]-preeop[8])*fixf;
        double ddeps   = preeop[9]+(nexteop[9]-preeop[9])*fixf;

        x_pole  = x_pole/Arcs;
        y_pole  = y_pole/Arcs;
        ddpsi   = ddpsi/Arcs;
        ddeps   = ddeps/Arcs;

        salida[0] = UT1_UTC;
        salida[1] = TAI_UTC;
        salida[2] = x_pole;
        salida[3] = y_pole;
        salida[4] = ddpsi;
        salida[5] = ddeps;

    }
    else if (interp == 'n')
    {
        double mj = (floor(Mjd_UTC));
        double nop = 20026;

        double aux[13];

        for (int i=0; i<nop; i++)
        {
            if (mj==eop[i][3])
            {
                for (int j=0; j<13; j++)
                {
                    aux[j] = eop[i][j] * 1.0;
                }

                break;
            }
        }

        double UT1_UTC = aux[6];
        double TAI_UTC = aux[12];
        double x_pole  = aux[4]/Arcs;
        double y_pole  = aux[5]/Arcs;
        double ddpsi   = aux[8]/Arcs;
        double ddeps   = aux[9]/Arcs;

        salida[0] = UT1_UTC;
        salida[1] = TAI_UTC;
        salida[2] = x_pole;
        salida[3] = y_pole;
        salida[4] = ddpsi;
        salida[5] = ddeps;
    }
}
