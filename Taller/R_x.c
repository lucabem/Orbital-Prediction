
#include "..\R_x.h"
#include <math.h>
#include "..\MatlabUtilities.h"

void R_x( long double angle,  long double matrizRotada[3][3])
{
     long double c = cos(angle);
     long double s = sin(angle);

    zeros(3,3, matrizRotada);

    matrizRotada[0][0] = 1.0;  matrizRotada[0][1] = 0.0;  matrizRotada[0][2] = 0.0;
    matrizRotada[1][0] = 0.0;  matrizRotada[1][1] = c;  matrizRotada[1][2] = s;
    matrizRotada[2][0] = 0.0;  matrizRotada[2][1] = -1.0*s;  matrizRotada[2][2] = c;
}
