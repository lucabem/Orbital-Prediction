
#include "PoleMatrix.h"
#include "R_x.h"
#include "R_y.h"
#include "MatlabUtilities.h"

void  PoleMatrix ( long double xp,  long double yp,  long double matrizTransformada[3][3])
{
     long double matrizTY[3][3], matrizTX[3][3];

    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            matrizTransformada[i][j] = 0.0;
            matrizTY[i][j] = 0.0;
            matrizTX[i][j] = 0.0;
        }
    }

    R_y(-xp, matrizTY);
    R_x(-yp, matrizTX);

    multiplicacion(3, 3, 3, 3, matrizTY, matrizTX, matrizTransformada);
}
