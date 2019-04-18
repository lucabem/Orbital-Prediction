
#include "PoleMatrix.h"
#include "R_x.h"
#include "R_y.h"
#include "MatlabUtilities.h"

void  PoleMatrix (double xp, double yp, double matrizTransformada[3][3])
{
     double matrizTY[3][3], matrizTX[3][3];

     R_y(xp, matrizTY); R_x(yp, matrizTX);

     multiplicacion(3, 3, 3, 3, matrizTY, matrizTX, matrizTransformada);
}
