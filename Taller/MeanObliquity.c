
#include "MeanObliquity.h"
#include "Constantes.h"

 double MeanObliquity( double Mjd_TT)
{

     double T = (Mjd_TT-MJD_J2000)/36525;

    return  RAD*(23.43929111-(46.8150+(0.00059-0.001813*T)*T)*T/3600);

}
