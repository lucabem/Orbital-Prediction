
#include "Anglesg.h"
#include "Constantes.h"
#include "MatlabUtilities.h"
#include "Gibbs.h"
#include "HGibbs.h"
#include "lambert_gooding.h"
#include "rv2coe.h"
#include "angl.h"
#include <stdio.h>
#include <string.h>

/**
    Accion que resuelve el problema de la determinacion de orbita usando tres
    vistas opticas.

    Inputs:
        Alpha1      - right ascension at t1    [rad]
        Alpha2      - right ascension at t2    [rad]
        Alpha3      - right ascension at t3    [rad]
        Delta1      - declination at t1         [rad]
        Delta2      - declination at t2         [rad]
        Delta3      - declination at t3         [rad]
        JD1         - Modified julian date of t1
        JD2         - Modified julian date of t2
        JD3         - Modified julian date of t3
        RS1         - ijk site1 position vector   [m]
        RS2         - ijk site2 position vector   [m]
        RS3         - ijk site3 position vector   [m]

    Outputs:
        vectorR2    - ijk position vector at t2   [m]
        V2          - ijk velocity vector at t2   [m/s]

*/
void anglesg( double Alpha1,  double Alpha2,  double Alpha3,  double Delta1,  double Delta2,  double Delta3,
              double JD1,  double JD2,  double JD3,
              double RS1[3],  double RS2[3],  double RS3[3], double vectorR2[3], double V2[3])
{
    double Mu = 398600.4418e9;
    double Rad = 180/PI;

    double RhoMat[3][1], CMat[3][1];

    double R1[3][1], R2[3][1], R3[3][1];
    zeros(3, 1, R1);
    zeros(3, 1, R2);
    zeros(3, 1, R3);

    double Tau1 = (JD1-JD2)*86400;
    double Tau3 = (JD3-JD2)*86400;

    double f1, g1, f3, g3;

    double L1[3], L2[3], L3[3];

    L1[0] = cos(Delta1)*cos(Alpha1);
    L1[1] = cos(Delta1)*sin(Alpha1);
    L1[2] = sin(Delta1);

    L2[0] = cos(Delta2)*cos(Alpha2);
    L2[1] = cos(Delta2)*sin(Alpha2);
    L2[2] = sin(Delta2);

    L3[0] = cos(Delta3)*cos(Alpha3);
    L3[1] = cos(Delta3)*sin(Alpha3);
    L3[2] = sin(Delta3);

    double LMatIi[3][3], RSMat[3][3];
    for (int i=0; i<3; i++)
    {
        LMatIi[i][0] = L1[i];
        LMatIi[i][1] = L2[i];
        LMatIi[i][2] = L3[i];
        RSMat[i][0] = RS1[i];
        RSMat[i][1] = RS2[i];
        RSMat[i][2] = RS3[i];
    }


    double D = det(3, LMatIi);

    double LMatI[3][3];

    LMatI[0][0] = ( L2[1]*L3[2]-L2[2]*L3[1])/D;
    LMatI[1][0] = (-L1[1]*L3[2]+L1[2]*L3[1])/D;
    LMatI[2][0] = ( L1[1]*L2[2]-L1[2]*L2[1])/D;
    LMatI[0][1] = (-L2[0]*L3[2]+L2[2]*L3[0])/D;
    LMatI[1][1] = ( L1[0]*L3[2]-L1[2]*L3[0])/D;
    LMatI[2][1] = (-L1[0]*L2[2]+L1[2]*L2[0])/D;
    LMatI[0][2] = ( L2[0]*L3[1]-L2[1]*L3[0])/D;
    LMatI[1][2] = (-L1[0]*L3[1]+L1[1]*L3[0])/D;
    LMatI[2][2] = ( L1[0]*L2[1]-L1[1]*L2[0])/D;

    double LIR[3][3];

    multiplicacion(3, 3, 3, 3, LMatI, RSMat, LIR);

    double a1  = Tau3/(Tau3 - Tau1);
    double a1u = (Tau3*((Tau3-Tau1)*(Tau3-Tau1) - Tau3*Tau3 ))/(6.0*(Tau3 - Tau1));
    double a3  = -Tau1 / (Tau3 - Tau1);
    double a3u = -(Tau1*((Tau3-Tau1)*(Tau3-Tau1) - Tau1*Tau1 ))/(6.0*(Tau3 - Tau1));

    double D1 = LIR[1][0]*a1 - LIR[1][1] + LIR[1][2]*a3;
    double D2 = LIR[1][0]*a1u + LIR[1][2]*a3u;

    double L2DotRS = dot(3, L2, RS2);
    double magRS2 = Norma(RS2);


    double Poly[16];
    Poly [ 0]=  1.0;  // r2^8th variable!!!!!!!!!!!!!!
    Poly [ 1]=  0.0;
    Poly [ 2]=  - (D1*D1 + 2.0*D1*L2DotRS + magRS2*magRS2);
    Poly [ 3]=  0.0;
    Poly [ 4]=  0.0;
    Poly [ 5]=  -2.0*Mu* (L2DotRS*D2 + D1*D2);
    Poly [ 6]=  0.0;
    Poly [ 7]=  0.0;
    Poly [ 8]=  -Mu*Mu*D2*D2;
    Poly [ 9]=  0.0;
    Poly [10]=  0.0;
    Poly [11]=  0.0;
    Poly [12]=  0.0;
    Poly [13]=  0.0;
    Poly [14]=  0.0;
    Poly [15]=  0.0;


    double zeror[15], zeroi[15];

    double resultado[15];

    for (int i=0; i<15; i++)
    {
        resultado[i] = 70;
    }


    int pos = raicesPolinomiales(15, Poly, zeror, zeroi, resultado);


    double BigR2 = 0.0;

    for (int i=0; i<15; i++)
    {
        if ((resultado[i]) > BigR2)
        {
            BigR2 = resultado[i];
        }
    }



    double u = Mu/(BigR2*BigR2*BigR2);

    double c1 = a1+a1u*u;
    double c3 = a3+a3u*u;
    CMat[0][0] = -c1;
    CMat[1][0] = 1.0;
    CMat[2][0] = -c3;

    multiplicacion(3, 3, 3, 1, LIR, CMat, RhoMat);

    double Rhoold2 = - RhoMat[1][0];

    double Rho2 = 999999e3;
    double ll = 0;


    //variables while

    double angulos[2], vecSalida[3], V1[3];
    double salidaRv2coe[11];
    double copa, theta, theta1;
    double p, a, ecc, incl, omega, argp, Nu, m, l, ArgPer;
    double magR2;

    double vectorR1[3] = {R1[0][0], R1[1][0], R1[2][0]};

    double vectorR3[3] = {R3[0][0], R3[1][0], R3[2][0]};


    while ( (fabs(Rhoold2-Rho2)>1e-12) && (ll<=2) )
    {

        ll = ll+1;
        Rho2 = Rhoold2;

        // reset now that inside while loop
        //---------- Now form the three position vectors ----------
        for (int i = 0; i<3; i++)
        {
            vectorR1[i] =  RhoMat[0][0]*L1[i]/c1 + RS1[i];
            vectorR2[i] = -RhoMat[1][0]*L2[i]    + RS2[i];
            vectorR3[i] =  RhoMat[2][0]*L3[i]/c3 + RS3[i];
        }



        char error[12] = "          ok";
        copa = gibbs(vectorR1, vectorR2, vectorR3, vecSalida, angulos, error);
        theta = angulos[0];
        theta1 = angulos[1];

        if ( strcmp(error,"ok") != 0 && (copa < 1/RAD) )
        {
            //--- HGibbs to get middle vector ----
            copa = hgibbs(vectorR1, vectorR2, vectorR3, JD1, JD2, JD3, vecSalida, angulos, error);
            theta = angulos[0];
            theta1 = angulos[1];
        }

        lambert_gooding(vectorR1, vectorR2, (JD2-JD1)*86400, Mu, false, 1, V1, V2);

        rv2coe( vectorR2,V2, salidaRv2coe);

        p = salidaRv2coe[0];
        a = salidaRv2coe[1];
        ecc = salidaRv2coe[2];
        incl = salidaRv2coe[3];
        omega = salidaRv2coe[4];
        argp = salidaRv2coe[5];
        Nu = salidaRv2coe[6];
        m = salidaRv2coe[7];
        u = salidaRv2coe[8];
        l = salidaRv2coe[9];
        ArgPer = salidaRv2coe[10];

        magR2 = Norma(vectorR2);

        if ( ll <= 2 )
        {
            //--- Now get an improved estimate of the f and g series --
            //       .or. can the analytic functions be found now??
            double U = Mu/(pow(magR2,3));
            double RDot = dot(3, vectorR2, V2)/magR2;
            double UDot = (-3.0*Mu*RDot)/(pow(magR2,4));

            double TauSqr= Tau1*Tau1;
            f1 =  1.0 - 0.5*U*TauSqr -(1.0/6.0)*UDot*TauSqr*Tau1
                  + (1.0/24.0) * U*U*TauSqr*TauSqr
                  + (1.0/30.0)*U*UDot*TauSqr*TauSqr*Tau1;
            g1 = Tau1 - (1.0/6.0)*U*Tau1*TauSqr - (1.0/12.0) *  UDot*TauSqr*TauSqr + (1.0/120.0)*U*U*TauSqr*TauSqr*Tau1 + (1.0/120.0)*U*UDot*TauSqr*TauSqr*TauSqr;
            TauSqr = Tau3*Tau3;
            f3 =  1.0 - 0.5*U*TauSqr -(1.0/6.0)*UDot*TauSqr*Tau3
                  + (1.0/24.0) * U*U*TauSqr*TauSqr
                  + (1.0/30.0)*U*UDot*TauSqr*TauSqr*Tau3;
            g3 = Tau3 - (1.0/6.0)*U*Tau3*TauSqr - (1.0/12.0) *
                 UDot*TauSqr*TauSqr
                 + (1.0/120.0)*U*U*TauSqr*TauSqr*Tau3
                 + (1.0/120.0)*U*UDot*TauSqr*TauSqr*TauSqr;
        }
        else
        {
            //-------- Now use exact method to find f and g -----------
            double Theta = angl(vectorR1, vectorR2);
            double Theta1 = angl(vectorR2, vectorR3);

            double magR1 = Norma(vectorR1);
            double magR3 = Norma(vectorR3);

            f1 = 1.0 - ( (magR1*(1.0 - cos(Theta))/p ) );
            g1 = ( magR1*magR2*sin(-theta) ) / sqrt(p);  // - ANGLE because backwards!!
            f3 = 1.0 - ( (magR3*(1.0 - cos(Theta1))/p ) );
            g3 = ( magR3*magR2*sin(theta1) )/sqrt(p);
        }

        c1 =  g3/(f1*g3 - f3*g1);
        c3 = -g1/(f1*g3 - f3*g1);




//----- Solve for all three ranges via matrix equation ----
        CMat[0][0] = -c1;
        CMat[1][0] = 1.0;
        CMat[2][0] = -c3;

        multiplicacion(3, 3, 3, 1, LIR, CMat, RhoMat);

        Rhoold2 =  RhoMat[1][0] * (-1);

    }


    for (int i = 0; i<3; i++)
    {
        vectorR1[i] =   RhoMat[0][0] * L1[i]/c1 + RS1[i];
        vectorR2[i] = - RhoMat[1][0] * L2[i] + RS2[i];
        vectorR3[i] =   RhoMat[2][0] * L3[i]/c3 + RS3[i];
    }



}


