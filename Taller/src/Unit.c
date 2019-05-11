
#include "..\MatlabUtilities.h"
#include "..\Unit.h"

void unit ( long double vect[3],  long double vectNormalizado[3])
{
     long double normaVect = Norma(vect);
     long double small = 0.000001;

    if ( normaVect > small)
        for (int i=0; i<3; i++)
            vectNormalizado[i] = vect[i]/normaVect;
    else
        for (int i=0; i<3; i++)
            vectNormalizado[i] = 0.0;
}

