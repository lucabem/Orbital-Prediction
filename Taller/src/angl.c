#include "..\MatlabUtilities.h"

/**
    Funcion que dado dos vectores de dimension 3, devuelve el angulo
    formado por los dos vectores

    Salida: theta [-pi, pi]
*/

 double angl( double vec1[],  double vec2[])
{
     double small     = 0.00000001;
     double undefined = 999999.1;

     double normaVec1 = Norma(vec1);
     double normaVec2 = Norma(vec2);
     double theta = undefined;

    if ( (normaVec1*normaVec2) > (small*small))
    {
         double temp = dot(3, vec1, vec2)/(normaVec1*normaVec2);
        if (fabs(temp) > 1)
            temp = Sing(temp);
        theta = acos(temp);
    }

    return theta;
}
