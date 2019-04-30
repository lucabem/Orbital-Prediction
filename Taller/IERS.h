#ifndef IERS_H_INCLUDED
#define IERS_H_INCLUDED

void IERS(long double(*eop)[13], long double Mjd_UTC, char interp, long double salida[6]);

#endif // IERS_H_INCLUDED
