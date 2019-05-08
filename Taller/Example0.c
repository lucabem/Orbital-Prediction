
#include "lambert_gooding.h"
#include <assert.h>
#include "MatlabUtilities.h"

int main()
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

}
