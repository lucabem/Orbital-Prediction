
#include "PrecMatrix.h"
#include "R_x.h"
#include "R_y.h"
#include "R_z.h"
#include "Constantes.h"
#include "MatlabUtilities.h"

void PrecMatrix(double Mjd_1, double Mjd_2, double matrizTransformada[3][3])
{
    double T  = (Mjd_1-MJD_J2000)/36525;
    double dT = (Mjd_2-Mjd_1)/36525;

    // Precession angles
    double zeta  =  ( (2306.2181+(1.39656-0.000139*T)*T)+ ((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/ARCS;
    double z     =  zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/ARCS;
    double theta =  ( (2004.3109-(0.85330+0.000217*T)*T)- ((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/ARCS;

    // Precession matrix
    double matrizX[3][3], matrizY[3][3], matrizZ[3][3];

    R_z(-z, matrizX) ; R_y(theta, matrizY) ; R_z(-zeta, matrizZ);

    double resultadoXY[3][3];
    multiplicacion(3, 3, 3, 3, matrizX, matrizY, resultadoXY);

    multiplicacion(3, 3, 3, 3, resultadoXY, matrizZ, matrizTransformada);

}
