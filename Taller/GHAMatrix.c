#include "GHAMatrix.h"
#include "R_z.h"
#include "Gast.h"

void ghaMatrix (double Mjd_UT1, double matrizTransformada[3][3], double (*eop)[13])
{
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            matrizTransformada[i][j] = 0.0;
        }
    }

    double g = gast(Mjd_UT1, eop);
    R_z(g, matrizTransformada);
}
