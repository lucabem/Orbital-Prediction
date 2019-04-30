#ifndef LAMBERT_GOODING_H_INCLUDED
#define LAMBERT_GOODING_H_INCLUDED

#include <stdbool.h>

void lambert_gooding( double r1, double r2,  double tof,  double mu, bool _way, int multi_revs,  double v1[3],  double v2[3]);

void vlamb( double gm, double r1[3],  double r2[3] ,  double th,  double tdelt,  double salida_n_vri_vti_vrf_vtf[5]);

void tlamb( double m, double q, double qsqfm1,  double x,  double n,  double salida_t_dt_d2t_d3t[4]);

 double d8rt( double x);

void xlamb ( double m,  double q,  double qsqfm1,  double tin,  double salida_n_x_xpl[3]);
#endif // LAMBERT_GOODING_H_INCLUDED
