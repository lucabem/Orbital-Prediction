
#include "MatlabUtilities.h"
#include "Unit.h"

void unit (double vect[3], double vectNormalizado[3])
{
    double normaVect = Norma(vect);
    double small = 0.000001;

    if ( normaVect > small)
        for (int i=0; i<3; i++)
            vectNormalizado[i] = vect[i]/normaVect;
    else
        for (int i=0; i<3; i++)
            vectNormalizado[i] = 0.0;
}

