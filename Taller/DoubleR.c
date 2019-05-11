#include "DoubleR.h"
#include <math.h>
#include "MatlabUtilities.h"
#include "Constantes.h"
#include <stdio.h>
#include <stdlib.h>

/**
    Accion que realiza el trabajo de iteracion para el algoritmo
    de long double-r
*/
void doubler ( long double cc1,  long double cc2,
               long double magrsite1,  long double magrsite2,
               long double magr1in,  long double magr2in,
               long double los1[3],  long double los2[3],  long double los3[3],
               long double rsite1[3],  long double rsite2[3],  long double rsite3[3],
               long double t1,  long double t3, char direct,
               //salida
               long double r2[3],  long double r3[3],  long double sal_f1_f2_q1_magr1_magr2_a_deltae32[7]
             )
{
    long double mu = 398600.4418e9;

    long double rho1 = (-cc1 + sqrt(pow(cc1,2)-4*(pow(magrsite1,2)-pow(magr1in,2))))/2;
    long double rho2 = (-cc2 + sqrt(pow(cc2,2)-4*(pow(magrsite2,2)-pow(magr2in,2))))/2;

    long double r1[3];

    r1[0] = rho1*los1[0] + rsite1[0];
    r1[1] = rho1*los1[1] + rsite1[1];
    r1[2] = rho1*los1[2] + rsite1[2];


    r2[0] = rho2*los2[0] + rsite2[0];
    r2[1] = rho2*los2[1] + rsite2[1];
    r2[2] = rho2*los2[2] + rsite2[2];

    long double magr1 = Norma(r1);
    long double magr2 = Norma(r2);

    long double w[3] = {0.0, 0.0, 0.0};
    if (direct == 'y')
    {

        cross(r1, r2, w);
        w[0] = w[0] / (magr1*magr2);
        w[1] = w[1] / (magr1*magr2);
        w[2] = w[2] / (magr1*magr2);


    }
    else
    {
        cross(r1, r2, w);
        w[0] = -w[0] / (magr1*magr2);
        w[1] = -w[1] / (magr1*magr2);
        w[2] = -w[2] / (magr1*magr2);
    }

    long double rho3 = -dot(3, rsite3, w) / dot(3, los3, w);

    r3[0] = rho3*los3[0] + rsite3[0];
    r3[1] = rho3*los3[1] + rsite3[1];
    r3[2] = rho3*los3[2] + rsite3[2];
    long double magr3 = Norma(r3);


    long double cosdv21 = dot(3, r2,r1)/(magr2*magr1);
    long double sindv21 = sqrt(1 - pow(cosdv21,2));
    long double dv21 = atan2(sindv21,cosdv21);

    long double cosdv31 = dot(3, r3,r1)/(magr3*magr1);
    long double sindv31 = sqrt(1 - pow(cosdv31,2));
    long double dv31 = atan2(sindv31,cosdv31);

    long double cosdv32 = dot(3, r3,r2)/(magr3*magr2);
    long double sindv32 = sqrt(1 - pow(cosdv32,2));

    long double c1, c3, p;

    if ( dv31 > PI)
    {
        c1 = (magr2*sindv32)/(magr1*sindv31);
        c3 = (magr2*sindv21)/(magr3*sindv31);
        p = (c1*magr1+c3*magr3-magr2)/(c1+c3-1);
    }
    else
    {
        c1 = (magr1*sindv31)/(magr2*sindv32);
        c3 = (magr1*sindv21)/(magr3*sindv32);
        p = (c3*magr3-c1*magr2+magr1)/(-c1+c3+1);
    }

    long double ecosv1 = p/magr1-1;
    long double ecosv2 = p/magr2-1;
    long double ecosv3 = p/magr3-1;

    long double esinv2;
    if ( ! fequal(dv21, PI) )
    {
        esinv2 = (-cosdv21*ecosv2+ecosv1)/sindv21;
    }
    else
    {
        esinv2 = (cosdv32*ecosv2-ecosv3)/sindv31;
    }

    long double e = sqrt(pow(ecosv2, 2) + pow(esinv2, 2));
    long double a = p / (1 - pow(e, 2));


    long double n, s, c, sinde32, cosde32, deltae32,
           sinde21, cosde21, deltae21, deltam32, deltam12;

    long double sindh32, sindh21, deltah32, deltah21;

    if (pow(e, 2) < 1)
    {
        n = sqrt(mu/(pow(a,3)));

        s = magr2/p*sqrt(1-pow(e,2))*esinv2;
        c = magr2/p*(pow(e,2)+ecosv2);

        sinde32 = magr3/sqrt(a*p)*sindv32-magr3/p*(1-cosdv32)*s;
        cosde32 = 1-magr2*magr3/(a*p)*(1-cosdv32);
        deltae32 = atan2(sinde32,cosde32);

        sinde21 = magr1/sqrt(a*p)*sindv21+magr1/p*(1-cosdv21)*s;
        cosde21 = 1-magr2*magr1/(a*p)*(1-cosdv21);
        deltae21 = atan2(sinde21,cosde21);

        deltam32 = deltae32+2*s*pow((sin(deltae32/2)),2)-c*sin(deltae32);
        deltam12 = -deltae21+2*s*pow((sin(deltae21/2)),2)+c*sin(deltae21);
    }
    else
    {
        n = sqrt(mu/(-pow(a,3)));

        s = magr2/p*sqrt(pow(e,2)-1)*esinv2;
        c = magr2/p*(pow(e,2)+ecosv2);

        sindh32 = magr3/sqrt(-a*p)*sindv32-magr3/p*(1-cosdv32)*s;
        sindh21 = magr1/sqrt(-a*p)*sindv21+magr1/p*(1-cosdv21)*s;

        deltah32 = log( sindh32 + sqrt(pow(sindh32,2) +1) );
        deltah21 = log( sindh21 + sqrt(pow(sindh21,2) +1) );

        deltam32 = -deltah32+2*s*pow((sinh(deltah32/2)),2)+c*sinh(deltah32);
        deltam12 = deltah21+2*s*pow((sinh(deltah21/2)),2)-c*sinh(deltah21);
        deltae32 = deltah32;
    }

    long double f1 = t1 - deltam12/n;
    long double f2 = t3 - deltam32/n;

    long double q1 = sqrt(pow(f1,2) + pow(f2,2));

    sal_f1_f2_q1_magr1_magr2_a_deltae32[0] = f1;
    sal_f1_f2_q1_magr1_magr2_a_deltae32[1] = f2;
    sal_f1_f2_q1_magr1_magr2_a_deltae32[2] = q1;
    sal_f1_f2_q1_magr1_magr2_a_deltae32[3] = magr1;
    sal_f1_f2_q1_magr1_magr2_a_deltae32[4] = magr2;
    sal_f1_f2_q1_magr1_magr2_a_deltae32[5] = a;
    sal_f1_f2_q1_magr1_magr2_a_deltae32[6] = deltae32;



}
