
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "MatlabUtilities.h"
#include "angl.h"
#include "lambert_gooding.h"
#include "Unit.h"
#include "timeDiff.h"
#include "R_x.h"
#include "R_y.h"
#include "R_z.h"
#include "PrecMatrix.h"
#include "Position.h"
#include "PoleMatrix.h"
#include "NutAngles.h"
#include "NewtonNu.h"
#include "Mjday.h"
#include "MeanObliquity.h"
#include "Frac.h"
#include "EqnEquinox.h"
#include "gmst.h"
#include "HGibbs.h"
#include "Gibbs.h"
#include "IERS.h"
#include "Gast.h"
#include "GHAMatrix.h"
#include "NutMatrix.h"
#include "DoubleR.h"
#include "AngleDr.h"
#include "Anglesg.h"

//Tiene que haber 27 Test_void
void Test_Angl();
void Test_Unit();
void Test_TimeDiff();
void Test_Rx();
void Test_Ry();
void Test_Rz();
void Test_PrecMatrix();
void Test_Position();
void Test_PoleMatrix();
void Test_NutAngles();
void Test_Newtonnu();
void Test_Mjday();
void Test_MeanObliquity();
void Test_Lambert();
void Test_Frac();
void Test_EqnEquinox();
void Test_Gmst();
void Test_HGibbs();
void Test_Gibbs();
void Test_IERS();
void Test_Gast();
void Test_GHAMatrix();
void Test_NutMatrix();
void Test_DoubleR();
void Test_AngleDr();
void Test_Anglesg();
void Test_Rv2coe();
void Test_Funciones();

void Test_Funciones()
{

    Test_Angl();
    Test_Unit();
    Test_TimeDiff();
    Test_Rx();
    Test_Ry();
    Test_Rz();
    Test_PrecMatrix();
    Test_Position();
    Test_PoleMatrix();
    Test_NutAngles();
    Test_Newtonnu();
    Test_Mjday();
    Test_MeanObliquity();
    Test_Lambert();
    Test_Frac();
    Test_EqnEquinox();
    Test_Gmst();
    Test_HGibbs();
    Test_Gibbs();
    Test_IERS();
    Test_Gast();
    Test_GHAMatrix();
    Test_NutMatrix();
    Test_DoubleR();
    Test_AngleDr();
    Test_Anglesg();

}

void Test_Angl()
{
    double r1[3] = {20387627.071752905846, 1865163.696333978558, -109943.68855587858707};
    double r2[3] = {20435422.352154411376, 1070699.4467182455119, 1012905.4914336511865};
    double theta = angl(r1, r2);

    assert(fequal(theta, 0.067208822931428643943) == true);

    double r3[3] = {20435422.352154411376, 1070699.4467182455119, 1012905.4914336511865};
    double r4[3] = {20398157.066625602543, 271778.61586978775449, 2131538.3954207645729};
    theta = angl(r3, r4);

    assert(fequal(theta, 0.067084233489783415272) == true);
}

void Test_Unit()
{
    double r1[3] = {-9341115904217.03125, 16158801718408.738281, 29720729511155.265625};
    double unit_r1[3];

    unit(r1, unit_r1);

    assert( (fequal(unit_r1[0],-0.26616375823995613858) == true));
    assert( (fequal(unit_r1[1], 0.46042543932939555829) == true));
    assert( (fequal(unit_r1[2], 0.84685610856739723662) == true));

}

void Test_TimeDiff()
{
    double diferencia[5];
    timeDiff(0.25802269087559637217, 34.0, diferencia);

    assert( fequal(diferencia[0], -33.741977309124401074) == true);
    assert( fequal(diferencia[1], -15.0) == true);
    assert( fequal(diferencia[2], -14.741977309124401074) == true);
    assert( fequal(diferencia[3], 66.184) == true);
    assert( fequal(diferencia[4], 15.0) == true);

    double diferencia2[5];
    timeDiff(0.25801238588253655459, 34.0, diferencia2);

    assert( fequal(diferencia2[0], -33.741987614117462613) == true);
    assert( fequal(diferencia2[1], -15.0) == true);
    assert( fequal(diferencia2[2], -14.741987614117462613) == true);
    assert( fequal(diferencia2[3], 66.184) == true);
    assert( fequal(diferencia2[4], 15.0) == true);
}

void Test_Rx()
{
    double m[3][3];
    R_x(-0.0000025678608802699049133, m);
    assert( fequal(m[0][0], 1.0) == true);
    assert( fequal(m[0][1], 0.0) == true);
    assert( fequal(m[0][2], 0.0) == true);
    assert( fequal(m[1][0], 0.0) == true);
    assert( fequal(m[1][1], 0.99999999999670308171) == true);
    assert( fequal(m[1][2], -0.0000025678608802670830231) == true);
    assert( fequal(m[2][0], 0.0) == true);
    assert( fequal(m[2][1], 0.0000025678608802670830231) == true);
    assert( fequal(m[2][2], 0.99999999999670308171) == true);


    R_x(0.40907147047232295112, m);
    assert( fequal(m[0][0], 1.0) == true);
    assert( fequal(m[0][1], 0.0) == true);
    assert( fequal(m[0][2], 0.0) == true);

    assert( fequal(m[1][0], 0.0) == true);
    assert( fequal(m[1][1], 0.91749054793879758485) == true);
    assert( fequal(m[1][2], 0.39775758250844805985) == true);

    assert( fequal(m[2][0], 0.0) == true);
    assert( fequal(m[2][1], -0.39775758250844805985) == true);
    assert( fequal(m[2][2], 0.91749054793879758485) == true);
}

void Test_Ry()
{
    double m[3][3];
    R_y(-0.000000075789200806792994845, m);
    assert( fequal(m[0][0], 1.0) == true);
    assert( fequal(m[0][1], 0.0) == true);
    assert( fequal(m[0][2], 0.000000075789200806792928671) == true);
    assert( fequal(m[1][0], 0.0) == true);
    assert( fequal(m[1][1], 1.0) == true);
    assert( fequal(m[1][2], 0.0) == true);
    assert( fequal(m[2][0], -0.000000075789200806792928671) == true);
    assert( fequal(m[2][1], 0.0) == true);
    assert( fequal(m[2][2], 1.0) == true);


    R_y(0.00091335104880993782557, m);
    assert( fequal(m[0][0], 0.99999958289495982644) == true);
    assert( fequal(m[0][1], 0.0) == true);
    assert( fequal(m[0][2], -0.00091335092182215902624) == true);

    assert( fequal(m[1][0], 0.0) == true);
    assert( fequal(m[1][1], 1.0) == true);
    assert( fequal(m[1][2], 0.0) == true);

    assert( fequal(m[2][0], 0.00091335092182215902624) == true);
    assert( fequal(m[2][1], 0.0) == true);
    assert( fequal(m[2][2], 0.99999958289495982644) == true);
}

void Test_Rz()
{
    double m[3][3];
    R_z(2.2594339459743442156, m);

    assert( fequal(m[0][0], -0.63548586103589599361) == true);
    assert( fequal(m[0][1], 0.7721125050298471848) == true);
    assert( fequal(m[0][2], 0.0) == true);

    assert( fequal(m[1][0], -0.7721125050298471848) == true);
    assert( fequal(m[1][1], -0.63548586103589599361) == true);
    assert( fequal(m[1][2], 0.0) == true);

    assert( fequal(m[2][0], 0.0) == true);
    assert( fequal(m[2][1], 0.0) == true);
    assert( fequal(m[2][2], 1.0) == true);
}

void Test_PrecMatrix()
{
    double pm[3][3];

    PrecMatrix(51544.5, 54977.68155855324585, pm);

    assert( fequal(pm[0][0], 0.99999737378108022323) == true);
    assert( fequal(pm[0][1], -0.0021019566973474981993) == true);
    assert( fequal(pm[0][2], -0.00091335041738156567178) == true);

    assert( fequal(pm[1][0], 0.0021019566973333333147) == true);
    assert( fequal(pm[1][1], 0.99999779088612039679) == true);
    assert( fequal(pm[1][2], -0.00000095992828239738784503) == true);

    assert( fequal(pm[2][0], 0.00091335041741416394482) == true);
    assert( fequal(pm[2][1], -0.0000009598972654118064989) == true);
    assert( fequal(pm[2][2], 0.99999958289495982644) == true);
}

void Test_Position()
{
    double pos[3];

    Position(-1.5047233973021472142, 0.53358904023671438477, 0.0, pos);
    assert(fabs(pos[0] - 362889.51475075335475) <  pow(10, -9));
    assert(fabs(pos[1] - (-5484262.3610134748742)) < pow(10, -9) );
    assert(fabs(pos[2]- 3225167.7284776144661)< pow(10, -9));
}

void Test_PoleMatrix()
{
    double pm[3][3];
    PoleMatrix(0.000000075789200806792994845021316663997, 0.0000025677719304258054046645378537539, pm);

    assert( fequal(pm[0][0], 1.0) == true);
    assert( fequal(pm[0][1], 0.00000000000019460938246087380948292199795477) == true);
    assert( fequal(pm[0][2], 0.000000075789200806543080420911895518643) == true);

    assert( fequal(pm[1][0], 0.0) == true);
    assert( fequal(pm[1][1], 0.99999999999670330375067806016887) == true);
    assert( fequal(pm[1][2], -0.0000025677719304229835144007601521743) == true);

    assert( fequal(pm[2][0], -0.000000075789200806792928670572312421783) == true);
    assert( fequal(pm[2][1], 0.0000025677719304229758911042348634712) == true);
    assert( fequal(pm[2][2], 0.99999999999670041717081403476186) == true);


    PoleMatrix(0.000000075972419978477384316826678551521, 0.0000025678608802699049133298671993009, pm);

    assert( fequal(pm[0][0], 1.0) == true);
    assert( fequal(pm[0][1], 0.00000000000019460938246087380948292199795477) == true);
    assert( fequal(pm[0][2], 0.000000075972419978226834618006816680913) == true);

    assert( fequal(pm[1][0], 0.0) == true);
    assert( fequal(pm[1][1], 0.99999999999670308170607313513756) == true);
    assert( fequal(pm[1][2], -0.0000025678608802670830230660894977213) == true);

    assert( fequal(pm[2][0], -0.000000075972419978477304907487873460864) == true);
    assert( fequal(pm[2][1], 0.0000025678608802670753997695642090182) == true);
    assert( fequal(pm[2][2], 0.99999999999670019512620910973055) == true);
}

void Test_NutAngles()
{
    double angulos[2];
    NutAngles(54977.667669664253481, angulos);
    assert( fequal(angulos[0], 0.000064869338739384949338) == true );
    assert( fequal(angulos[1], 0.000022305133670726570955) == true );

    NutAngles(54977.68155855324585, angulos);
    assert( fequal(angulos[0], 0.000064882163989950450574) == true );
    assert( fequal(angulos[1], 0.000022305027665978605734) == true );

}

void Test_Newtonnu()
{
    double newt[2];
    newtonnu(0.082533106173374254366, 0.1812003110698841013, newt);
    assert( fequal(newt[0], 0.16688393775668777796) == true);
    assert( fequal(newt[1], 0.1531743313690979158) == true);

    newtonnu(0.079129107778551310837, 0.19026723747469370673, newt);
    assert( fequal(newt[0], 0.17584036930451643621) == true);
    assert( fequal(newt[1], 0.16199787056861986168) == true);

}

void Test_Mjday()
{
    assert( fabs (Mjday(2011, 1, 4, 13, 0, 46.5) -  55565.54220486106351) < 0.0000000001);
    assert( fabs(Mjday(2011, 1, 4, 13, 10, 46.5)- 55565.549149305559695) < 0.0000000001);

}

void Test_MeanObliquity()
{
    assert( fequal(MeanObliquity(55565.542970879585482), 0.40906781750982068591 ) == true);
    assert( fequal(MeanObliquity(55565.546443102066405), 0.4090678174882443896  ) == true);
}

void Test_Lambert()
{
    double salida_tlamb[4];
    tlamb(0.0, 0.804611564466231, 0.352600230327204, 0.0, 0.0, salida_tlamb);
    assert(fequal(salida_tlamb[0], 2.227109718426920) == true );
    assert(fequal(salida_tlamb[1], 0.000000000000000) == true );
    assert(fequal(salida_tlamb[2], 0.000000000000000) == true );
    assert(fequal(salida_tlamb[3], 0.000000000000000) == true );

    assert(fequal(d8rt(0.00554491061386597),0.522380333950713) == true);
    assert(fequal(d8rt(0.137320935252477),0.780220027563195) == true);

    double salida_xlamb[3];
    xlamb(0.0, 0.804611564466231, 0.352600230327204, 0.913438160983312, salida_xlamb);
    assert(fequal(salida_xlamb[0], 1.0) == true);
    assert(fequal(salida_xlamb[1], 0.62563831524107) == true);
    assert(fequal(salida_xlamb[2], 0.0) == true);

    xlamb(1.0, 0.804611564466231, 0.352600230327204, 0.913438160983312, salida_xlamb);
    assert(fequal(salida_xlamb[0], 0.0) == true);
    assert(fequal(salida_xlamb[1], 0.0) == true);
    assert(fequal(salida_xlamb[2], 0.0) == true);


    double vri[2][1], vti[2][1], vrf[2][1], vtf[2][1];
    int n;

    n = vlamb(398600441800000, 41206854.7330225, 40716998.0363724, 0.012666404089535, 299.999982118607, vri, vti, vrf, vtf);
    assert(n == 1);

    assert(fabs(vri[0][0]- (-1469.867334575)) < pow(10, -9));
    assert(fabs(vti[0][0]- (1577.84211894393)) < pow(10, -9));
    assert(fabs(vrf[0][0]- (-1527.41527291996)) < pow(10, -9));
    assert(fabs(vtf[0][0]- (1596.82476907768)) < pow(10, -9));

    assert(fequal(vri[1][0], 0.0) == true);
    assert(fequal(vti[1][0], 0.0) == true);
    assert(fequal(vrf[1][0], 0.0) == true);
    assert(fequal(vtf[1][0], 0.0) == true);

    n = vlamb(398600441800000, 41206854.7330225, 40716998.0363724, 6.29585171126912, 299.999982118607, vri, vti, vrf, vtf);
    assert(fequal(vri[0][0], 0.0) == true);
    assert(fequal(vti[0][0], 0.0) == true);
    assert(fequal(vrf[0][0], 0.0) == true);
    assert(fequal(vtf[0][0], 0.0) == true);

    assert(fequal(vri[0][0], 0.0) == true);
    assert(fequal(vti[0][0], 0.0) == true);
    assert(fequal(vrf[0][0], 0.0) == true);
    assert(fequal(vtf[0][0], 0.0) == true);

    double lg_v1[3];
    double lg_v2[3];

    double r1[3] = {8794276.5809840224683, 404708.19494348927401, 2543973.8056371584535};
    double r2[3] = {8330586.9962050709873, 3762923.0821287548169, 572416.99634550604969};

    lambert_gooding(r1, r2, 600.00000447034835815, 398600441800000, false, 1, lg_v1, lg_v2);
    assert( fabs(lg_v1[0]- 591.41567971816994032) < pow(10, -5));
    assert( fabs(lg_v1[1]- 5838.8636504518080983)< pow(10, -5));
    assert( fabs(lg_v1[2]- (-2988.6398832611826037))< pow(10, -5));

    assert( fabs(lg_v2[0] - (-2113.6537704040365497))< pow(10, -4));
    assert( fabs(lg_v2[1] - 5180.3929960643363302)< pow(10, -4));
    assert( fabs(lg_v2[2] - (-3480.8307130797988975)) < pow(10, -4));
}

void Test_Frac()
{
    assert( fequal(Frac(4.5), 0.5) == true);
    assert( fequal(Frac(0.4444), 0.4444) == true);
    assert( fequal(Frac(-4.6666), 0.3334) == true);


}

void Test_EqnEquinox()
{
    assert( fequal(eqnEquinox(2323232), -7.61848099380581e-05) == true);
    assert( fequal(eqnEquinox(0), 2.55262176011142e-05) == true);
    assert( fequal(eqnEquinox(-5555550), -6.97398129880609e-05) == true);
}

void Test_Gmst()
{
    assert( fequal(gmst(5550), 2.20092031214172) == true);
    assert( fequal(gmst(0), 0.973208148169487) == true);
    assert( fequal(gmst(-6700), 5.09502762284038) == true);

}

void Test_HGibbs()
{
    double r1[3] = {20387627.0717529, 1865163.69633398, -109943.688555879};
    double r2[3] = {20435422.3521544, 1070699.44671825, 1012905.49143365};
    double r3[3] = {20398157.0666256, 271778.615869788, 2131538.39542076};
    double vectVel[3], angulos[2];
    char error[12] = "";

    double copa = hgibbs(r1, r2, r3, 55565.9044073611, 55565.9078795835, 55565.9113518056, vectVel, angulos, error );

    assert (fabs( vectVel[0] - 17.4309174129376) < pow(10,-5));
    assert (fabs( vectVel[1] + 2657.49386520154) < pow(10,-5));
    assert (fabs( vectVel[2] - 3738.39266893706) < pow(10,-5));

    assert (fabs( angulos[0] - 0.0672088229314286) < pow(10,-12));
    assert (fabs( angulos[1] - 0.0670842334897834) < pow(10,-12));

    assert( fequal(copa, -7.87130777224476e-16) == true);
    assert( strcmp(error, "   angl > 1ø") == 0);

}

void Test_Gibbs()
{
    double r1[3] = {20387627.0717529, 1865163.69633398, -109943.688555879};
    double r2[3] = {20435422.3521544, 1070699.44671825, 1012905.49143365};
    double r3[3] = {20398157.0666256, 271778.615869788, 2131538.39542076};
    double vectVel[3] = {0.0, 0.0, 0.0};
    double angulos[2] = {0.0, 0.0};
    char error[12] = "";

    double copa = gibbs(r1, r2, r3, vectVel, angulos, error);

    assert (fabs( vectVel[0] - 17.4448460090308) < pow(10,-5));
    assert (fabs( vectVel[1] + 2659.68695020331) < pow(10,-5));
    assert (fabs( vectVel[2] -  3741.47770465728) < pow(10,-5));

    assert (fabs( angulos[0] - 0.0672088229314286) < pow(10,-12));
    assert (fabs( angulos[1] - 0.0670842334897834) < pow(10,-12));

    assert( fequal(copa, -7.87130777224476e-16) == true);
}

void Test_IERS()
{
    double (*eop)[13] = malloc(sizeof( double[20026][13]));

    FILE* fid = fopen("eop19620101.txt","rt");

    int v1, v2, v3, v4, v13;
    float v5, v6, v7, v8, v9, v10, v11, v12;
    if (fid == NULL)
    {
        exit(EXIT_FAILURE);
    }

    int fila = 0;
    while( fscanf(fid,"%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12, &v13) != EOF)
    {
        eop[fila][0] = v1;
        eop[fila][1] = v2;
        eop[fila][2] = v3;
        eop[fila][3] = v4;
        eop[fila][4] = v5;
        eop[fila][5] = v6;
        eop[fila][6] = v7;
        eop[fila][7] = v8;
        eop[fila][8] = v9;
        eop[fila][9] = v10;
        eop[fila][10] = v11;
        eop[fila][11] = v12;
        eop[fila][12] = v13;

        fila++;
    }
    fclose(fid);

    double salida[6];

    IERS(eop, 55565.5422048611, 'l', salida);

    assert( fabs(salida[0]+0.141248008109364) < pow(10, -6));
    assert( fabs(salida[1]-34) < pow(10, -6));
    assert( fabs(salida[2]- 5.78543831547586e-07) < pow(10, -5));
    assert( fabs(salida[3]- 9.72285830548474e-07)< pow(10, -5));
    assert( fabs(salida[4]+3.23115555150206e-07)< pow(10, -5));
    assert( fabs(salida[5]+3.0657888529631e-08)< pow(10, -5));
}

void Test_Gast()
{
    double (*eop)[13] = malloc(sizeof( double[20026][13]));

    FILE* fid = fopen("eop19620101.txt","rt");

    int v1, v2, v3, v4, v13;
    float v5, v6, v7, v8, v9, v10, v11, v12;
    if (fid == NULL)
    {
        exit(EXIT_FAILURE);
    }

    int fila = 0;
    while( fscanf(fid,"%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12, &v13) != EOF)
    {
        eop[fila][0] = v1;
        eop[fila][1] = v2;
        eop[fila][2] = v3;
        eop[fila][3] = v4;
        eop[fila][4] = v5;
        eop[fila][5] = v6;
        eop[fila][6] = v7;
        eop[fila][7] = v8;
        eop[fila][8] = v9;
        eop[fila][9] = v10;
        eop[fila][10] = v11;
        eop[fila][11] = v12;
        eop[fila][12] = v13;

        fila++;
    }
    fclose(fid);


    assert( fabs(gast(55565.5422032263, eop) - 5.21832550097282 ) < pow(10, -9));
    assert( fabs(gast(55565.5491476707, eop) - 5.26207819997555 ) < pow(10, -9));

}

void Test_GHAMatrix()
{
    double (*eop)[13] = malloc(sizeof( double[20026][13]));

    FILE* fid = fopen("eop19620101.txt","rt");

    int v1, v2, v3, v4, v13;
    float v5, v6, v7, v8, v9, v10, v11, v12;
    if (fid == NULL)
    {
        exit(EXIT_FAILURE);
    }

    int fila = 0;
    while( fscanf(fid,"%d %d %d %d %f  %f  %f  %f  %f  %f  %f  %f   %d", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12, &v13) != EOF)
    {
        eop[fila][0] = v1;
        eop[fila][1] = v2;
        eop[fila][2] = v3;
        eop[fila][3] = v4;
        eop[fila][4] = v5;
        eop[fila][5] = v6;
        eop[fila][6] = v7;
        eop[fila][7] = v8;
        eop[fila][8] = v9;
        eop[fila][9] = v10;
        eop[fila][10] = v11;
        eop[fila][11] = v12;
        eop[fila][12] = v13;

        fila++;
    }
    fclose(fid);

    double m[3][3];

    ghaMatrix(55565.5422032263, m, eop);

    assert( fabs(m[0][0] - 0.484626846950978)  < pow(10, -9));
    assert( fabs(m[0][1] + 0.874720995068915 ) < pow(10, -9));
    assert( fabs(m[0][2] - 0.0 )< pow(10, -9));

    assert( fabs(m[1][0] - 0.874720995068915 )< pow(10, -9));
    assert( fabs(m[1][1] - 0.484626846950978 )< pow(10, -9));
    assert( fabs(m[1][2] - 0.0 )< pow(10, -9));

    assert( fabs(m[2][0] - 0.0 )< pow(10, -9));
    assert( fabs(m[2][1] - 0.0 )< pow(10, -9));
    assert( fabs(m[2][2] - 1.0 )< pow(10, -9));

}

void Test_NutMatrix()
{
    double m[3][3];

    NutMatrix(55565.5499153241, m);

    assert( fabs(m[0][0] - 0.999999996195242)  < pow(10, -9));
    assert( fabs(m[0][1] + 8.00351553835159e-05) < pow(10, -9));
    assert( fabs(m[0][2] +3.46971108657198e-05 )< pow(10, -9));

    assert( fabs(m[1][0] - 8.00351835718741e-05  )< pow(10, -9));
    assert( fabs(m[1][1] - 0.999999996796856 )< pow(10, -9));
    assert( fabs(m[1][2] - 8.11024524267579e-07 )< pow(10, -9));

    assert( fabs(m[2][0] - 3.46970458441061e-05  )< pow(10, -9));
    assert( fabs(m[2][1] + 8.1380151086029e-07 )< pow(10, -9));
    assert( fabs(m[2][2] - 0.999999999397726 )< pow(10, -9));
}

void Test_DoubleR()
{
    double los1[3] = {0.148851929355429, -0.840681876398112, -0.520669843396865};
    double los3[3] = {0.173629378644301, -0.832793490773584, -0.525649922093347};
    double los2[3] = {0.161099913661194, -0.836869774686652, -0.523159438445173};


    double rsite1[3] = {-4625314.68025011, -2963495.20415729, 3230277.01672135};
    double rsite2[3] = {-4559382.005115, -3064041.00297889, 3230203.98499169};
    double rsite3[3] = {-4491265.71743158, -3163120.47843576, 3230128.5447305};

    double r2[3], r3[3], sal_f1_f2_q1_magr1_magr2_a_deltae32[7];


    doubler(241923.732943492, 279551.107139144, 6372639.11744252, 6372639.11744252, 12820055.37, 13457869.07, los1, los2, los3, rsite1, rsite2, rsite3, -300.000022351742, 299.999982118607, 'y',r2, r3, sal_f1_f2_q1_magr1_magr2_a_deltae32);

    assert( fabs( r2[0]+2672181.15108884) < pow(10, -7));
    assert( fabs( r2[1]+12867530.76039)< pow(10, -7));
    assert( fabs( r2[2]+2898333.99239024)< pow(10, -7));

    assert( fabs( r3[0]+2326788.44033551)< pow(10, -7));
    assert( fabs( r3[1]+13544788.611832)< pow(10, -7));
    assert( fabs( r3[2]+3322664.0815058)< pow(10, -7));

    assert( fabs( sal_f1_f2_q1_magr1_magr2_a_deltae32[0]+192.733796845969)< pow(10, -7));
    assert( fabs( sal_f1_f2_q1_magr1_magr2_a_deltae32[1]-184.377089731178)< pow(10, -7));
    assert( fabs( sal_f1_f2_q1_magr1_magr2_a_deltae32[2]-266.723129226549)< pow(10, -7));
    assert( fabs( sal_f1_f2_q1_magr1_magr2_a_deltae32[3]-12820055.37)< pow(10, -7));
    assert( fabs( sal_f1_f2_q1_magr1_magr2_a_deltae32[4]-13457869.07)< pow(10, -7));
    assert( fabs( sal_f1_f2_q1_magr1_magr2_a_deltae32[5]-370142443.928873)< 0.5);
    assert( fabs( sal_f1_f2_q1_magr1_magr2_a_deltae32[6]-0.00869694889817696)< pow(10, -7));

}

void Test_AngleDr()
{
    double rsite1[3] = {4950990.3382646, 256563.116260381, 3999465.34658133};
    double rsite2[3] = {4935037.85913036, 472703.320202615, 3999475.70573182};
    double rsite3[3] = {4909646.95198536, 687938.936915757, 3999494.94894739};

    double r2[3] = {0.0, 0.0, 0.0};
    double v2[3] = {0.0, 0.0, 0.0};

    anglesdr(5.39901096780381, 6.26556833768239, 0.732191050658823, 0.0115360853036144,
             -0.360600766331881, -0.640322905313176, 54977.6669036457, 54977.6738480902, 54977.6807925347, rsite1, rsite2, rsite3, r2, v2 );

    assert( fabs((r2[0] - 8794276.58098403))< 0.01);
    assert( fabs(r2[1] - 404.708194943489*1000)< 0.01);
    assert( fabs(r2[2] - 2543.97380563715*1000) < 0.01);

    assert( fabs(  1000*(v2[0] - 0.591415680079019*1000)< 10));
    assert( fabs(  1000*(v2[1] - 5.83886365198696*1000))< 10);
    assert( fabs(  1000*(v2[2] + 2.98863988398241*1000)) < 10);

}

void Test_Anglesg()
{
    double rs1[3] = {5270137.35006701, -1572248.25164427, 3219350.41084204};
    double rs2[3] = {5303269.31336066, -1456667.74823777, 3219314.05463487};
    double rs3[3] = {5333865.06903306, -1340390.16746883, 3219280.50120838};
    double r2[3], v2[3];

    anglesg(0.223578422509726, 0.165492119674102, 0.106613437358074, -0.211533905341713, -0.142837745983216, -0.0716736910623991, 55565.9044073611, 55565.9078795835, 55565.9113518056, rs1, rs2, rs3, r2, v2);

    assert( fabs( r2[0]/1000 - 20486.511511687) < pow(10,-2));
    assert( fabs( r2[1]/1000 - 1079.23234124354) < pow(10,-2));
    assert( fabs( r2[2]/1000 - 1005.45621711092) < pow(10,-2));

    assert( fabs( v2[0]/1000 - 0.0168797950290674) < pow(10,-2));
    assert( fabs( v2[1]/1000 + 2.65408002932635) < pow(10,-2));
    assert( fabs( v2[2]/1000 - 3.73412004615382) < pow(10,-2));

}
