
#include "lambert_gooding.h"
#include "MatlabUtilities.h"
#include <math.h>
#include "Constantes.h"
#include "Unit.h"
#include <stdio.h>
#include <stdlib.h>
void lambert_gooding( long double r1[3], long double r2[3],  long double tof,  long double mu, bool long_way, int multi_revs,  long double v1[3],  long double v2[3])
{
    long double r1mag = Norma(r1);
    long double r2mag = Norma(r2);

    int sols = 0;

    if ( r1mag==0.0 || r2mag==0.0 || mu<=0.0 || tof<=0.0 )
    {
        printf("Error in solve_lambert_gooding: invalid input\n");
        return ;
    }

    int solution_exists[5] = {20, 20, 20, 20, 20};
    long double dr       = r1mag - r2mag;
    long double r1r2     = r1mag*r2mag;
    long double r1hat[3] = {r1[0]/r1mag, r1[1]/r1mag, r1[2]/r1mag};
    long double r2hat[3] = {r2[0]/r2mag, r2[1]/r2mag, r2[2]/r2mag};
    long double r1xr2[3];
    long double all_vt1[3][2];
    long double all_vt2[3][2];

    cross(r1,r2, r1xr2);

    if (r1xr2[0] == 0.0 && r1xr2[1] == 0.0 && r1xr2[2] == 0.0)
    {
        r1xr2[0] = 0.0 ;
        r1xr2[1] = 0.0 ;
        r1xr2[2] = 1.0;
    }

    long double r1xr2_hat[3];
    unit(r1xr2,r1xr2_hat);

    long double pa = acos(fmax(-1.0, fmin(1.0,dot(3,r1hat,r2hat))));

    for (int i=0; i<=multi_revs; i++)
    {
        int num_revs = i;
        long double ta, rho[3];

        if (long_way)
        {
            ta    =  num_revs * 2*PI + (2*PI - pa);
            rho[0]   = -r1xr2_hat[0];
            rho[1]   = -r1xr2_hat[1];
            rho[2]   = -r1xr2_hat[2];

        }
        else
        {
            ta    = num_revs * 2*PI + pa;
            rho[0]   = r1xr2_hat[0];
            rho[1]   = r1xr2_hat[1];
            rho[2]   = r1xr2_hat[2];
        }

        long double etai[3], etaf[3];
        cross(rho, r1hat, etai);
        cross(rho, r2hat, etaf);

        long double vri[2][1], vti[2][1], vrf[2][1], vtf[2][1];
        int n = vlamb(mu,r1mag,r2mag,ta,tof, vri, vti, vrf, vtf);

        long double vt1[3][2];
        long double vt2[3][2];

        switch(n)
        {
        case 1:
            vt1[0][0] = vri[0][0]*r1hat[0]+ vti[0][0]*etai[0];
            vt1[1][0] = vri[0][0]*r1hat[1]+ vti[0][0]*etai[1];
            vt1[2][0] = vri[0][0]*r1hat[2]+ vti[0][0]*etai[2];

            vt2[0][0] = vrf[0][0]*r2hat[0]+ vtf[0][0]*etaf[0];
            vt2[1][0] = vrf[0][0]*r2hat[1]+ vtf[0][0]*etaf[1];
            vt2[2][0] = vrf[0][0]*r2hat[2]+ vtf[0][0]*etaf[2];

        case 2:
            vt1[0][0] = vri[0][0]*r1hat[0]+ vti[0][0]*etai[0];
            vt1[1][0] = vri[0][0]*r1hat[1]+ vti[0][0]*etai[1];
            vt1[2][0] = vri[0][0]*r1hat[2]+ vti[0][0]*etai[2];

            vt2[0][0] = vrf[0][0]*r2hat[0]+ vtf[0][0]*etaf[0];
            vt2[1][0] = vrf[0][0]*r2hat[1]+ vtf[0][0]*etaf[1];
            vt2[2][0] = vrf[0][0]*r2hat[2]+ vtf[0][0]*etaf[2];

            vt1[0][1] = vri[0][0]*r1hat[0]+ vti[0][0]*etai[0];
            vt1[1][1] = vri[0][0]*r1hat[1]+ vti[0][0]*etai[1];
            vt1[2][1] = vri[0][0]*r1hat[2]+ vti[0][0]*etai[2];

            vt2[0][1] = vrf[0][0]*r2hat[0]+ vtf[1][0]*etaf[0];
            vt2[1][1] = vrf[0][0]*r2hat[1]+ vtf[1][0]*etaf[1];
            vt2[2][1] = vrf[0][0]*r2hat[2]+ vtf[1][0]*etaf[2];
        }

        if ( i==0 && n==1)
        {
            all_vt1[0][0] = vt1[0][0];
            all_vt1[1][0] = vt1[1][0];
            all_vt1[2][0] = vt1[2][0];

            all_vt2[0][0] = vt2[0][0];
            all_vt2[1][0] = vt2[1][0];
            all_vt2[2][0] = vt2[2][0];

            sols++;
            solution_exists[0] = 1;
        }
        else
        {

            switch(n)
            {
            case 1:
                all_vt1[0][2*i-1] = vt1[0][0];
                all_vt1[1][2*i-1] = vt1[1][0];
                all_vt1[2][2*i-1] = vt1[2][0];

                all_vt2[0][2*i-1] = vt2[0][0];
                all_vt2[1][2*i-1] = vt2[1][0];
                all_vt2[2][2*i-1] = vt2[2][0];

                sols++;
                solution_exists[2*i-1] = 1;

            case 2:
                all_vt1[0][2*i] = vt1[0][0];
                all_vt1[1][2*i] = vt1[1][0];
                all_vt1[2][2*i] = vt1[2][0];

                all_vt2[0][2*i] = vt2[0][0];
                all_vt2[1][2*i] = vt2[1][0];
                all_vt2[2][2*i] = vt2[2][0];

                sols++;
                solution_exists[2*i-1] = 1;

                all_vt1[0][2*i] = vt1[0][0];
                all_vt1[1][2*i] = vt1[1][0];
                all_vt1[2][2*i] = vt1[2][0];

                all_vt2[0][2*i] = vt2[0][0];
                all_vt2[1][2*i] = vt2[1][0];
                all_vt2[2][2*i] = vt2[2][0];

                sols++;
                solution_exists[2*i] = 1;
            }
        }
    }

    int n_solutions = sols;


    int k = 0;

    for (int i=0; i<5; i++)
    {
        if ( fequal(solution_exists[i], 1) )
        {
            v1[0] = all_vt1[0][i];
            v1[1] = all_vt1[1][i];
            v1[2] = all_vt1[2][i];

            v2[0] = all_vt2[0][i];
            v2[1] = all_vt2[1][i];
            v2[2] = all_vt2[2][i];

            k = k + 1;

        }
    }


}


int vlamb( long double gm, long double r1,  long double r2,  long double th,  long double tdelt, long double vri[2][1], long double vti[2][1], long double vrf[2][1], long double vtf[2][1] )
{
    zeros(2, 1, vri);
    zeros(2, 1, vti);
    zeros(2, 1, vrf);
    zeros(2, 1, vtf);

    long double thr2 = th;

    int m = 0;

    while (thr2 > 2*PI)
    {
        thr2 = thr2 - 2*PI;
        m = m + 1;

    }


    thr2 = thr2/2;

    long double r1mag = r1;
    long double r2mag = r2;
    long double dr = r1mag-r2mag;
    long double r1r2 = r1mag*r2mag;
    long double r1r2th = 4.0*r1r2*sin(thr2)*sin(thr2);
    long double csq    = dr*dr + r1r2th;
    long double c      = sqrt(csq);
    long double s      = (r1 + r2 + c)/2.0;
    long double gms    = sqrt(gm*s/2.0);
    long double qsqfm1 = c/s;
    long double q      = sqrt(r1r2)*cos(thr2)/s;

    long double rho, sig;
    if ( !(fequal(c,0.0)) )
    {
        rho = dr/c;
        sig = r1r2th/csq;

    }
    else
    {
        rho = 0.0;
        sig = 1.0;
    }


    long double t = 4.0*gms*tdelt/(s*s);

    long double salida[3];


    xlamb(m,q,qsqfm1,t, salida);

    long double n  = salida[0];
    long double x1 = salida[1];
    long double x2 = salida[2];


    long double x;

    for (int i=1; i<=n; i++)
    {

        if (i==1)
        {
            x = x1;
        }
        else
        {
            x = x2;
        }


        long double sal[4];

        tlamb(m, q, qsqfm1, x, -1, sal);
        long double qzminx = sal[1];
        long double qzplx  = sal[2];
        long double zplqx  = sal[3];

        long double vt2 = gms*zplqx*sqrt(sig);
        long double vr1 = gms*(qzminx - qzplx*rho)/r1;
        long double vt1 = vt2/r1;
        long double vr2 = -gms*(qzminx + qzplx*rho)/r2;
        vt2 = vt2/r2;


        vri[i-1][0] = vr1;
        vti[i-1][0] = vt1;
        vrf[i-1][0] = vr2;
        vtf[i-1][0] = vt2;
    }


    return n;
}

void tlamb( long double m, long double q, long double qsqfm1,  long double x,  long double n,  long double salida_t_dt_d2t_d3t[4])
{
    //printf("Entrada \n");
    //printf("\t m=%.15f q=%.15f qsqfm1=%.15f x=%.15f n=%.15f \n", m, q, qsqfm1, x, n);
    long double dt, d2t, d3t;
    long double sw = 0.4;
    long double t = 0;

    bool lm1 = (n==-1);
    bool l1 = (n>=1);
    bool l2 = (n>=2);
    bool l3 = (n==3);
    long double qsq = q*q;
    long double xsq = x*x;
    long double u = (1.0 - x)*(1.0 + x);

    if (!lm1)
    {
        // (needed if series, and otherwise useful when z = 0)
        dt = 0.0;
        d2t = 0.0;
        d3t = 0.0;
    }

    if (lm1 || m>0 || x<0.0 || fabs(u)>sw)
    {
        // direct computation (not series)
        long double y = sqrt(fabs(u));
        long double z = sqrt(qsqfm1 + qsq*xsq);
        long double qx = q*x;

        long double a, b, aa, bb, g, f, fg1sq;
        long double fg1,term, twoi1, told, qz, qz2;
        if (qx<=0.0)
        {
            a = z - qx;
            b = q*z - x;
        }

        if (qx<0.0 && lm1)
        {
            aa = qsqfm1/a;
            bb = qsqfm1*(qsq*u - xsq)/b;
        }
        if ((qx==0.0&&lm1) || (qx>0.0))
        {
            aa = z + qx;
            bb = q*z + x;
        }
        if (qx>0.0)
        {
            a = qsqfm1/aa;
            b = qsqfm1*(qsq*u - xsq)/bb;
        }
        if (!lm1)
        {
            if (qx*u>=0.0)
                g = x*z + q*u;
            else
                g = (xsq - qsq*u)/(x*z - q*u);
            f = a*y;

            if (x<=1.0)
            {
                t = m*PI + atan2(f, g);
            }
            else
            {
                if (f>sw)
                    t = log(f + g);
                else
                {
                    fg1 = f/(g + 1.0);
                    term = 2.0*fg1;
                    fg1sq = fg1*fg1;
                    t = term;
                    twoi1 = 1.0;
                    while(1)
                    {
                        twoi1 = twoi1 + 2.0;
                        term = term*fg1sq;
                        told = t;
                        t = t + term/twoi1;
                        if ( !fequal(t,told))
                        {
                            break;
                        }
                    } // (continue looping for inverse tanh)
                }
            }


            t = 2.0*(t/y + b)/u;

            if (l1 && !fequal(z,0.0))
            {
                qz = q/z;
                qz2 = qz*qz;
                qz = qz*qz2;
                dt = (3.0*x*t - 4.0*(a + qx*qsqfm1)/z)/u;
                if (l2)
                {
                    d2t = (3.0*t + 5.0*x*dt + 4.0*qz*qsqfm1)/u;
                }
                if (l3)
                {
                    d3t = (8.0*dt + 7.0*x*d2t - 12.0*qz*qz2*x*qsqfm1)/u;
                }
            }
        }
        else
        {
            dt = b;
            d2t = bb;
            d3t = aa;
        }
    }
    else
    {
        // compute by series
        long double u0i = 1.0;
        long double u1i, u2i, u3i;
        if (l1)
        {
            u1i = 1.0;
        }
        if (l2)
        {
            u2i = 1.0;
        }
        if (l3)
        {
            u3i = 1.0;
        }
        long double term = 4.0;
        long double tq = q*qsqfm1;
        int i = 0;

        long double tqsum, ttmold;

        if (q<0.5)
        {
            tqsum = 1.0 - q*qsq;
        }
        if (q>=0.5)
        {
            tqsum = (1.0/(1.0 + q) + q)*qsqfm1;
        }

        ttmold = term/3.0;
        t = ttmold*tqsum;

        while(1)
        {
            i = i + 1;
            int p = i;
            u0i = u0i*u;
            if (l1 && i>1)
            {
                u1i = u1i*u;
            }
            if (l2 && i>2)
            {
                u2i = u2i*u;
            }
            if (l3 && i>3)
            {
                u3i = u3i*u;
            }
            term = term*(p - 0.5)/p;
            tq = tq*qsq;
            tqsum = tqsum + tq;
            long double told = t;
            long double tterm = term/(2.0*p + 3.0);
            long double tqterm = tterm*tqsum;
            t = t - u0i*((1.5*p + 0.25)*tqterm/(p*p - 0.25)-ttmold*tq);
            ttmold = tterm;
            tqterm = tqterm*p;
            if (l1)
            {
                dt = dt + tqterm*u1i;
            }
            if (l2)
            {
                d2t = d2t + tqterm*u2i*(p - 1.0);
            }
            if (l3)
            {
                d3t = d3t + tqterm*u3i*(p - 1.0)*(p - 2.0);
            }
            if (i<n || !fequal(t,told))
            {
                // cycle
                break;
            }
        }
        if (l3)
        {
            d3t = 8.0*x*(1.5*d2t - xsq*d3t);
        }
        if (l2)
        {
            d2t = 2.0*(2.0*xsq*d2t - dt);
        }
        if (l1)
        {
            dt = -2.0*x*dt;
        }
        t = t/xsq;
    }

    salida_t_dt_d2t_d3t[0] = t;
    salida_t_dt_d2t_d3t[1] = dt;
    salida_t_dt_d2t_d3t[2] = d2t;
    salida_t_dt_d2t_d3t[3] = d3t;


//    printf("Salida \n");
//    printf("\t t=%.15f dt=%.15f d2t=%.15f d3t=%.15f\n", t, dt, d2t, d3t);
}



long double d8rt( long double x)
{
    return pow(x, 0.125);
}

void xlamb ( long double m,  long double q,  long double qsqfm1,  long double tin,  long double salida_n_x_xpl[3])
{

    long double tol = 3e-7;
    long double c0  = 1.7;
    long double c1  = 0.5;
    long double c2  = 0.03;
    long double c3  = 0.15;
    long double c41 = 1.0;
    long double c42 = 0.24;

    long double thr2 = atan2(qsqfm1, 2.0*q)/PI;

    long double xpl = 0;
    long double x = 0;
    int n = 0;

    long double salida_tmin_dt_d2t_d3t[4];

    long double t, t0, dt, d2t, d3t, tdiff, w, xm, tmin, tdiffm, d2t2, tdiff0;

    if ( m == 0)
    {

        n = 1;

        tlamb(m, q, qsqfm1, 0, 0, salida_tmin_dt_d2t_d3t);

        t0 = salida_tmin_dt_d2t_d3t[0];
        dt = salida_tmin_dt_d2t_d3t[1];
        d2t = salida_tmin_dt_d2t_d3t[2];
        d3t = salida_tmin_dt_d2t_d3t[3];

        tdiff = tin - t0;

        if ( tdiff <= 0.0)
        {
            x = t0*tdiff/(-4.0*tin);
        }
        else
        {
            x = -tdiff/(tdiff + 4.0);
            w = x + c0*sqrt(2.0*(1.0 - thr2));
            if (w<0.0)
                x = x - sqrt(d8rt(-w))*(x + sqrt(tdiff/(tdiff + 1.5*t0)));

            w = 4.0/(4.0 + tdiff);
            x = x*(1.0 + x*(c1*w - c2*x*sqrt(w)));
        }
    }
    else
    {

        xm = 1.0/(1.5*(m + 0.5)*PI);
        if (thr2<0.5)
            xm = d8rt(2.0*thr2)*xm;

        if (thr2>0.5)
            xm = (2.0 - d8rt(2.0 - 2.0*thr2))*xm;

        int i;
        for (i=1; i<=12; i++)
        {
            tlamb(m,q,qsqfm1,xm,3, salida_tmin_dt_d2t_d3t );

            tmin = salida_tmin_dt_d2t_d3t[0];
            dt   = salida_tmin_dt_d2t_d3t[1];
            d2t  = salida_tmin_dt_d2t_d3t[2];
            d3t  = salida_tmin_dt_d2t_d3t[3];

            if (d2t == 0.0)
            {
                break;
            }

            long double xmold = xm;
            xm = xm - dt*d2t/(d2t*d2t - dt*d3t/2.0);
            long double xtest = fabs(xmold/xm - 1.0);
            if (xtest<=tol)
                break;
        }

        if ( i > 12)
        {
            n = -1;
            salida_n_x_xpl[0] = n*1.0;
            salida_n_x_xpl[1] = x;
            salida_n_x_xpl[2] = xpl;
            return;
        }

        tdiffm = tin - tmin;

        if ( tdiffm < 0.0)
        {
            n = 0;
            salida_n_x_xpl[0] = n*1.0;
            salida_n_x_xpl[1] = x;
            salida_n_x_xpl[2] = xpl;
            return;

        }
        else if (tdiffm == 0.0)
        {
            x = xm;
            n = -1;
            salida_n_x_xpl[0] = n*1.0;
            salida_n_x_xpl[1] = x;
            salida_n_x_xpl[2] = xpl;
            return;
        }
        else
        {
            n = 3;

            if (d2t==0.0)
            {
                d2t = 6.0*m*PI;
            }


            x = sqrt(tdiffm/(d2t/2.0 + tdiffm/pow((1.0 - xm),2)));
            w = xm + x;
            w = w*4.0/(4.0 + tdiffm) + pow((1.0 - w),2);
            x = x*(1.0 - (1.0 + m + c41*(thr2 - 0.5))/(1.0 + c3*m)*x*(c1*w + c2*x*sqrt(w))) + xm;
            d2t2 = d2t/2.0;

            if (x>=1.0)
            {
                n = 1;
            }


            tlamb(m,q,qsqfm1,0.0,0, salida_tmin_dt_d2t_d3t);
            t0 = salida_tmin_dt_d2t_d3t[0];
            dt = salida_tmin_dt_d2t_d3t[1];
            d2t = salida_tmin_dt_d2t_d3t[2];
            d3t = salida_tmin_dt_d2t_d3t[3];

            tdiff0 = t0 - tmin;
            tdiff = tin - t0;

            if (tdiff<=0)
            {
                x = xm - sqrt(tdiffm/(d2t2 - tdiffm*(d2t2/tdiff0 - 1.0/pow(xm,2))));
            }
            else
            {
                x = -tdiff/(tdiff + 4.0);
                long double ij = 200;
                w = x + c0*sqrt(2.0*(1.0 - thr2));
                if (w<0.0)
                {
                    x = x - sqrt(d8rt(-w))*(x + sqrt(tdiff/(tdiff+1.5*t0)));
                }

                w = 4.0/(4.0 + tdiff);
                x = x*(1.0 + (1.0 + m + c42*(thr2 - 0.5))/(1.0 + c3*m)*x*(c1*w - c2*x*sqrt(w)));
                if (x<=-1.0)
                {
                    n = n - 1;
                    if (fequal(n,1))
                    {
                        x = xpl;
                    }

                }

            }
        }
    }

    for (int i=1; i<=3; i++)
    {
        tlamb(m,q,qsqfm1,x,2, salida_tmin_dt_d2t_d3t);
        t = salida_tmin_dt_d2t_d3t[0];
        dt = salida_tmin_dt_d2t_d3t[1];
        d2t = salida_tmin_dt_d2t_d3t[2];
        t = tin - t;

        if ( !fequal(dt,0.0) )
            x = x + t*dt/(dt*dt + t*d2t/2.0);
    }


    if (n!=3)
    {
        salida_n_x_xpl[0] = n*1.0;
        salida_n_x_xpl[1] = x;
        salida_n_x_xpl[2] = xpl;
        return;
    }
    else
    {
        n = 2;
        xpl = x;
        tlamb(m,q,qsqfm1,0.0,0, salida_tmin_dt_d2t_d3t);

        long double tdiff0 = salida_tmin_dt_d2t_d3t[0] - tmin;
        tdiff = tin - salida_tmin_dt_d2t_d3t[0];

        if (tdiff<=0)
            x = xm - sqrt(tdiffm/(d2t2 - tdiffm*(d2t2/tdiff0 - 1.0/pow(xm,2))));
        else
        {
            x = -tdiff/(tdiff + 4.0);
            long double ij = 200;
            w = x + c0*sqrt(2.0*(1.0 - thr2));
            if (w<0.0)
                x = x - sqrt(d8rt(-w))*(x + sqrt(tdiff/(tdiff+1.5*t0)));

            w = 4.0/(4.0 + tdiff);
            x = x*(1.0 + (1.0 + m + c42*(thr2 - 0.5))/(1.0 + c3*m)*x*(c1*w - c2*x*sqrt(w)));

            if ( x <= -1.0)
            {
                n = n - 1;
                if ( fequal(n, 1) )
                {
                    x = xpl;
                }
            }

        }

    }


    salida_n_x_xpl[0] = n*1.0;
    salida_n_x_xpl[1] = x;
    salida_n_x_xpl[2] = xpl;
}
