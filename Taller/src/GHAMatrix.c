#include "..\GHAMatrix.h"
#include "..\R_z.h"
#include "..\Gast.h"

/**
    Accion que transforma desde el verdadero equador y equinoccio al
    sistema de equador terraqueo y Greenwich meridiano.


    Entrada:
        long double Mjd_UT1:                   Fecha Juliana modificada
        long double matrizTransormada[3][3]:   Guarda el resultado de la transformacion
        long double (*eop)[13]:                Variable con datos almacenados
*/

void ghaMatrix ( long double Mjd_UT1,  long double matrizTransformada[3][3],  long double (*eop)[13])
{
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            matrizTransformada[i][j] = 0.0;
        }
    }

     long double g = gast(Mjd_UT1, eop);
    R_z(g, matrizTransformada);
}
