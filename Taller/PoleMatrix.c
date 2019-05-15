
#include "PoleMatrix.h"
#include "R_x.h"
#include "R_y.h"
#include "MatlabUtilities.h"

/**

     PoleMatrix: Transformación de pseudo-Tierra-fija a Tierra-fija
     de coordenadas para una fecha dada

     De entrada:
        Coordinador polar (xp, yp)

     De salida:
        PoleMat matriz de polo

    Last modified:   2015/08/12   M. Mahooti

*/
void  PoleMatrix ( double xp,  double yp,  double matrizTransformada[3][3])
{
     double matrizTY[3][3], matrizTX[3][3];

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
