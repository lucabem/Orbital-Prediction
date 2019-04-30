#include "Gibbs.h"
#include "Constantes.h"
#include "MatlabUtilities.h"
#include "Unit.h"
#include <math.h>
#include "angl.h"

long double gibbs (long double r1[3],
              long double r2[3],
              long double r3[3],
              long double vectVelocidadSal[3],
              long double angulos[2],
              char *error
             )
{
    long double small = 0.00000001;
    long double theta = 0;
    long double theta1 = 0;

    long double v2[3][1];
    zeros(3, 1, v2);
    error = "ok";

    long double normavect1 = Norma(r1);
    long double normavect2 = Norma(r2);
    long double normavect3 = Norma(r3);

    long double p[3], q[3], w[3];
    cross(r2, r3, p);
    cross(r3, r1, q);
    cross(r1, r2, w);

    long double pn[3], r1n[3];

    unit(p, pn);
    unit(r1, r1n);
    long double copa = asin(dot(3, pn, r1n));

    if(fabs(dot(3, r1n, pn)) > 0.017452406 )
    {
        error = "not coplanar";
    }

    long double vectorSuma[3] = {p[0]+q[0]+w[0], p[1]+q[1]+w[1], p[2]+q[2]+w[2] };
    long double normaSuma = Norma(vectorSuma);

    long double n[3] = {normavect1*p[0]+normavect2*q[0]+normavect3*w[0],
                   normavect1*p[1]+normavect2*q[1]+normavect3*w[1],
                   normavect1*p[2]+normavect2*q[2]+normavect3*w[2]
                  };

    long double normaN = Norma(n);
    long double nn[3], dn[3];

    unit(n, nn);
    unit(vectorSuma, dn);

    if ( ( fabs(normaSuma)<small ) || ( fabs(normaN)<small ) ||  ( dot(3, nn,dn) < small ) ){
         error= "impossible";
        vectVelocidadSal[0] = 0.0;
        vectVelocidadSal[1] = 0.0;
        vectVelocidadSal[2] = 0.0;
    }
    else
    {
        theta  = angl(r1,r2);
        theta1 = angl(r2,r3);


//      ----------- perform gibbs method to find v2 -----------
        long double r1mr2 = normavect1 - normavect2;
        long double r3mr1 = normavect3 - normavect1;
        long double r2mr3 = normavect2 - normavect3;

        long double s[3]  = {r1mr2*r3[0] + r3mr1*r2[0] + r2mr3*r1[0], r1mr2*r3[1] + r3mr1*r2[1] + r2mr3*r1[1], r1mr2*r3[2] + r3mr1*r2[2] + r2mr3*r1[2]};

        long double b[3];
        cross(vectorSuma,r2, b);

        long double l  = sqrt(GM_Earth/(normaSuma*normaN));
        long double tover2 = l/normavect2;

        long double v2[3] = {tover2 * b[0] + l * s[0], tover2 * b[1] + l * s[1], tover2 * b[2] + l * s[2]};

        vectVelocidadSal[0] = v2[0];
        vectVelocidadSal[1] = v2[1];
        vectVelocidadSal[2] = v2[2];
    }

    angulos[0] = theta;
    angulos[1] = theta1;

    return copa;
}
