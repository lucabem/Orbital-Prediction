#include "Anglesg.h"
#include "Constantes.h"
#include "MatlabUtilities.h"

void anglesg(double Alpha1, double Alpha2, double Alpha3, double Delta1, double Delta2, double Delta3,
             double JD1[3], double JD2[3], double JD3[3],
             double RS1[3], double RS2[3], double RS3[3])
{
    double Mu = 398600.4418e9;
    double Rad = 180/PI;

    double R1[3][1], R2[3][1], R3[3][1];
    zeros(3, 1, R1);
    zeros(3, 1, R2);
    zeros(3, 1, R3);

    double Tau1 = (JD1-JD2)*86400;
    double Tau3 = (JD3-JD2)*86400;


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
    double *resultado = raicesPolinomiales(Poly, 15, zeror, zeroi);

    double BigR2 = 0.0;

    for (int i=0; i<15; i++)
    {
        if ((resultado[i]) > BigR2)
        {
            BigR2 = resultado[i];
        }
    }
}
