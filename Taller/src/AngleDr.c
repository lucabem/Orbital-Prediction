#include "..\AngleDr.h"
#include "..\Constantes.h"
#include <math.h>
#include "..\MatlabUtilities.h"
#include "..\DoubleR.h"
#include "..\lambert_gooding.h"

/**
    Accion que resuelve el problema de la determinacion de orbita usando tres
    vistas opticas.

    Inputs:
        rtasc1      - right ascension at t1    [rad]
        rtasc2      - right ascension at t2    [rad]
        rtasc3      - right ascension at t3    [rad]
        decl1       - declination at t1         [rad]
        decl2       - declination at t2         [rad]
        decl3       - declination at t3         [rad]
        Mjd1        - Modified julian date of t1
        Mjd2        - Modified julian date of t2
        Mjd3        - Modified julian date of t3
        rsite1      - ijk site1 position vector   [m]
        rsite2      - ijk site2 position vector   [m]
        rsite3      - ijk site3 position vector   [m]

    Outputs:
        r2        - ijk position vector at t2   [m]
        v2        - ijk velocity vector at t2   [m/s]


*/
void anglesdr(double rtasc1, double rtasc2, double rtasc3,
              double decl1, double decl2, double decl3,
              double Mjd1, double Mjd2, double Mjd3,
              double rsite1[3], double rsite2[3], double rsite3[3],
              double r2[3], double v2[3])

{
    double magr1in = 2.01 * R_Earth;
    double magr2in = 2.11 * R_Earth;
    char direct  = 'y';


    double tol    = 1e-8*R_Earth;
    double pctchg = 5e-6;

    double t1 = (Mjd1 - Mjd2)*86400;
    double t3 = (Mjd3 - Mjd2)*86400;

    double los1[3] = {cos(decl1)*cos(rtasc1), cos(decl1)*sin(rtasc1), sin(decl1)};
    double los2[3] = {cos(decl2)*cos(rtasc2), cos(decl2)*sin(rtasc2), sin(decl2)};
    double los3[3] = {cos(decl3)*cos(rtasc3), cos(decl3)*sin(rtasc3), sin(decl3)};

    double magr1old  = 99999e3;
    double magr2old  = 99999e3;
    double magrsite1 = Norma(rsite1);
    double magrsite2 = Norma(rsite2);
    double magrsite3 = Norma(rsite3);


    double cc1 = 2*dot(3,los1,rsite1);
    double cc2 = 2*dot(3,los2,rsite2);


    int ll = 0;

    double r3[3], sal_f1_f2_q1_magr1_magr2_a_deltae32[7];


    while  ((fabs(magr1in-magr1old) > tol && fabs(magr2in-magr2old) > tol && ll<=3))
    {
        ll++;

        doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct, r2, r3, sal_f1_f2_q1_magr1_magr2_a_deltae32);

        double f1,f2,q1,magr1,magr2,a,deltae32;

        f1 = sal_f1_f2_q1_magr1_magr2_a_deltae32[0];
        f2 = sal_f1_f2_q1_magr1_magr2_a_deltae32[1];
        q1 = sal_f1_f2_q1_magr1_magr2_a_deltae32[2];
        magr1 = sal_f1_f2_q1_magr1_magr2_a_deltae32[3];
        magr2 = sal_f1_f2_q1_magr1_magr2_a_deltae32[4];
        a = sal_f1_f2_q1_magr1_magr2_a_deltae32[5];
        deltae32 = sal_f1_f2_q1_magr1_magr2_a_deltae32[6];


        //ouble v3[3];
        //lambert_gooding(r2, r3, (Mjd3-Mjd2)*86400, GM_Earth, false, 1, v2, v3);

        double f  = 1 - (a/magr2)*(1-cos(deltae32));
        double g  = t3 - sqrt(pow(a,3)/GM_Earth)*(deltae32-sin(deltae32));
        double v2[3];
        v2[0] = (r3[0]- f*r2[0])/g;
        v2[1] = (r3[1] - f*r2[1])/g;
        v2[2] = (r3[2] - f*r2[2])/g;




        double magr1o = magr1in;
        magr1in = (1+pctchg)*magr1in;
        double deltar1 = pctchg*magr1in;


        doubler(cc1,cc2,
                magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,
                rsite3,t1,t3,direct, r2, r3, sal_f1_f2_q1_magr1_magr2_a_deltae32);

        double f1delr1 = sal_f1_f2_q1_magr1_magr2_a_deltae32[0];
        double f2delr1 = sal_f1_f2_q1_magr1_magr2_a_deltae32[1];
        double q2 = sal_f1_f2_q1_magr1_magr2_a_deltae32[2];
        magr1 = sal_f1_f2_q1_magr1_magr2_a_deltae32[3];
        magr2 = sal_f1_f2_q1_magr1_magr2_a_deltae32[4];
        a = sal_f1_f2_q1_magr1_magr2_a_deltae32[5];
        deltae32 = sal_f1_f2_q1_magr1_magr2_a_deltae32[6];

        double pf1pr1 = (f1delr1-f1)/deltar1;
        double pf2pr1 = (f2delr1-f2)/deltar1;

        magr1in = magr1o;
        deltar1 = pctchg*magr1in;
        double magr2o = magr2in;
        magr2in = (1+pctchg)*magr2in;
        double deltar2 = pctchg*magr2in;



        doubler(cc1,cc2,
                magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,
                rsite3,t1,t3,direct, r2, r3, sal_f1_f2_q1_magr1_magr2_a_deltae32);

        double f1delr2 = sal_f1_f2_q1_magr1_magr2_a_deltae32[0];
        double f2delr2 = sal_f1_f2_q1_magr1_magr2_a_deltae32[1];
        double q3 = sal_f1_f2_q1_magr1_magr2_a_deltae32[2];
        magr1 = sal_f1_f2_q1_magr1_magr2_a_deltae32[3];
        magr2 = sal_f1_f2_q1_magr1_magr2_a_deltae32[4];
        a = sal_f1_f2_q1_magr1_magr2_a_deltae32[5];
        deltae32 = sal_f1_f2_q1_magr1_magr2_a_deltae32[6];


        double pf1pr2 = (f1delr2-f1)/deltar2;
        double pf2pr2 = (f2delr2-f2)/deltar2;



        magr2in = magr2o;
        deltar2 = pctchg*magr2in;


        double delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
        double delta1 = pf2pr2*f1 - pf1pr2*f2;
        double delta2 = pf1pr1*f2 - pf2pr1*f1;



        deltar1 = -delta1/delta;
        deltar2 = -delta2/delta;



        magr1old = magr1in;
        magr2old = magr2in;

/*
        if (fabs(deltar1) > magr1in*pctchg)
        {
            deltar1 = Sing(deltar1)*magr1in*pctchg;
        }

        if (fabs(deltar2) > magr2in*pctchg)
        {
            deltar2 = Sing(deltar2)*magr2in*pctchg;
        }
*/
        magr1in = magr1in + deltar1;
        magr2in = magr2in + deltar2;



    }

    doubler(cc1,cc2,magrsite1,
            magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct, r2, r3, sal_f1_f2_q1_magr1_magr2_a_deltae32);


    double f1,f2,q1,magr1,magr2,a,deltae32;

    f1 = sal_f1_f2_q1_magr1_magr2_a_deltae32[0];
    f2 = sal_f1_f2_q1_magr1_magr2_a_deltae32[1];
    q1 = sal_f1_f2_q1_magr1_magr2_a_deltae32[2];
    magr1 = sal_f1_f2_q1_magr1_magr2_a_deltae32[3];
    magr2 = sal_f1_f2_q1_magr1_magr2_a_deltae32[4];
    a = sal_f1_f2_q1_magr1_magr2_a_deltae32[5];
    deltae32 = sal_f1_f2_q1_magr1_magr2_a_deltae32[6];


    double v3[3];
    lambert_gooding(r2, r3, (Mjd3-Mjd2)*86400,GM_Earth, false, 1, v2, v3);

    double f  = 1 - a/magr2*(1-cos(deltae32));
    double g  = t3 - sqrt(pow(a,3)/GM_Earth)*(deltae32-sin(deltae32));


    v2[0] = (r3[0]- f*r2[0])/g;
    v2[1] = (r3[1] - f*r2[1])/g;
    v2[2] = (r3[2] - f*r2[2])/g;



}
