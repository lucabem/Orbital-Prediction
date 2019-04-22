

#include "NutMatrix.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"
#include "MatlabUtilities.h"

void NutMatrix(double Mjd_TT, double matrizNutacion[3][3])
{
    double ep = MeanObliquity(Mjd_TT);

    // Nutation in longitude and obliquity
    double angulosNut[2] = {0.0, 0.0};
    NutAngles(Mjd_TT, angulosNut);

    // Transformation from mean to true equator and equinox
    double A[3][3], B[3][3], C[3][3];
    R_x(-ep-angulosNut[1], A);
    R_z(-angulosNut[0], B);
    R_x(ep, C);

    double productoAB[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    multiplicacion(3, 3, 3, 3, A, B, productoAB );


    multiplicacion(3, 3, 3, 3, productoAB, C, matrizNutacion );

}
