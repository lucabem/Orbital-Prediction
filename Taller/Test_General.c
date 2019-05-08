
#include <assert.h>
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

void Test_Angl();
void Test_Unit();
void Test_TimeDiff();
void Test_Rx();
void Test_Ry();
void Test_Rz();
void Test_PrecMatrix();
void Test_Position();
void Test_PoleMatrix();

void Test_AngleDr();
void Test_DoubleR();
void Test_EqnEquinox();
void Test_Frac();
void Test_Gast();
void Test_GHAMatrix();
void Test_Gibbs();
void Test_Gmst();
void Test_HGibbs();
void Test_IERS();
void Test_Lambert();
void Test_MeanObliquity();
void Test_Mjday();
void Test_Newtonnu();
void Test_NutAngles();
void Test_NutMatrix();

int main()
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


