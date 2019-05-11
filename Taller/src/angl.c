#include "..\MatlabUtilities.h"

/**
    Funcion que dado dos vectores de dimension 3, devuelve el angulo
    formado por los dos vectores

    Salida: theta [-pi, pi]
*/

 long double angl( long double vec1[],  long double vec2[])
{
     long double small     = 0.00000001;
     long double undefined = 999999.1;

     long double normaVec1 = Norma(vec1);
     long double normaVec2 = Norma(vec2);
     long double theta = undefined;

    if ( (normaVec1*normaVec2) > (small*small))
    {
         long double temp = dot(3, vec1, vec2)/(normaVec1*normaVec2);
        if (fabs(temp) > 1)
            temp = Sing(temp);
        theta = acos(temp);
    }

    return theta;
}
