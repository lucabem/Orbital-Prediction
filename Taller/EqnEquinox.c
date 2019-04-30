#include "EqnEquinox.h"
#include "NutAngles.h"
#include "MeanObliquity.h"
#include "math.h"

long double eqnEquinox (long double fechaJulianaModificada )
{
    long double vectAngulos[2] = {0.0, 0.0};
    NutAngles(fechaJulianaModificada, vectAngulos);
    return vectAngulos[0] * cos(MeanObliquity(fechaJulianaModificada));
}
