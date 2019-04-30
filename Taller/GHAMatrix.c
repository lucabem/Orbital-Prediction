#include "GHAMatrix.h"
#include "R_z.h"
#include "Gast.h"

void ghaMatrix (long double Mjd_UT1, long double matrizTransformada[3][3], long double (*eop)[13])
{
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            matrizTransformada[i][j] = 0.0;
        }
    }

    long double g = gast(Mjd_UT1, eop);
    R_z(g, matrizTransformada);
}
