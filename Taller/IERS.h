#ifndef IERS_H_INCLUDED
#define IERS_H_INCLUDED

void IERS(double(*eop)[13], double Mjd_UTC, char interp, double salida[6]);

#endif // IERS_H_INCLUDED
