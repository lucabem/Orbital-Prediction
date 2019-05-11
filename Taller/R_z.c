
#include "R_z.h"
#include "MatlabUtilities.h"

void R_z( double angle,  double matrizRotada[3][3])
{
     double c = cos(angle);
     double s = sin(angle);

    zeros(3,3, matrizRotada);

    matrizRotada[0][0] =      c;  matrizRotada[0][1] =   s;  matrizRotada[0][2] = 0.0;
    matrizRotada[1][0] = -1.0*s;  matrizRotada[1][1] =   c;  matrizRotada[1][2] = 0.0;
    matrizRotada[2][0] =    0.0;  matrizRotada[2][1] = 0.0;  matrizRotada[2][2] = 1.0;
}
