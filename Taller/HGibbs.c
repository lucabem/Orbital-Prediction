#include "HGibbs.h"
#include "Constantes.h"
#include "MatlabUtilities.h"
#include <math.h>
#include "Unit.h"
#include "angl.h"

double hgibbs(double r1[3], double r2[3], double r3[3],
              double MJD1,double MJD2,double MJD3,
              double vectVel[3], double angulos[2], char* error)
{
    error = "ok";
    angulos[0] = 0.0, angulos[1] = 0.0;
    double magr1 = Norma(r1);
    double magr2 = Norma(r2);
    double magr3 = Norma(r3);

    for (int i=0; i<3; i++)
    {
        vectVel[i] = 0.0;
    }

    double tolAngle = 0.01745329251994;
    double dt21= (MJD2-MJD1)*86400;
    double dt31= (MJD3-MJD1)*86400;
    double dt32= (MJD3-MJD2)*86400;

    double p[3], pn[3], r1n[3];
    cross(r2, r3, p);
    unit(p, pn);
    unit(r1, r1n);
    double copa = asin(dot(3, pn, r1n));

    if ( fabs(dot(3, r1n, pn)) > 0.017452406 )
        error= "not coplanar";

    angulos[0] = angl(r1, r2);
    angulos[1] = angl(r2, r3);

    if ( (angulos[0] > tolAngle ) || (angulos[1] > tolAngle))
        error = "angle > 10";


    double term1= -dt32*( 1/(dt21*dt31) + GM_Earth/(12*magr1*magr1*magr1) );
    double term2= (dt32-dt21)*( 1/(dt21*dt32) + GM_Earth/(12*magr2*magr2*magr2) );
    double term3=  dt21*( 1/(dt32*dt31) + GM_Earth/(12*magr3*magr3*magr3) );

    vectVel[0] = term1*r1[0] + term2*r2[0] + term3*r3[0];
    vectVel[1] = term1*r1[1] + term2*r2[1] + term3*r3[1];
    vectVel[2] = term1*r1[2] + term2*r2[2] + term3*r3[2];

    return copa;
}
