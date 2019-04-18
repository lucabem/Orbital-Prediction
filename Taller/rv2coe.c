
#include "rv2coe.h"
#include "MatlabUtilities.h"

void rv2coe (double vectorPosicion[3], double vectorVelocidad[3], double vectorResultado[11])
{
    double mu = 398600.4418e9;
    double small = 1e-10;
    double undefined = 999999.1;

    double normaVectPosicion = Norma(vectorPosicion);
    double normaVectVelocidad = Norma(vectorVelocidad);

    double* hbar = cross(vectorPosicion, vectorVelocidad);
    double vHbar[3];
    for (int i=0; i<3; i++)
    {
        vHbar[i] = hbar[i];
    }

    double normaHbar = Norma(vHbar);

    if ( normaHbar > small)
    {
        double nbar[3];
        nbar[0] = -vHbar[1];
        nbar[1] =  vHbar[0];
        nbar[2] =   0.0;

        double normaNbar = Norma(nbar);
        double c1 = normaVectVelocidad*normaVectVelocidad - mu / normaVectPosicion;

        double rdotv = dot(3, vectorPosicion, vectorVelocidad);
        double ebar[3];
        for (int i=0; i<3; i++)
               ebar[i] = (c1*vectorPosicion[i] - rdotv*vectorVelocidad[i])/mu;

        double ecc = Norma( ebar );
    }
}
