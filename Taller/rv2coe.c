
#include "rv2coe.h"
#include "MatlabUtilities.h"
#include "Constantes.h"
#include <math.h>
#include "angl.h"
#include <string.h>
#include "NewtonNu.h"

void rv2coe (long double vectorPosicion[3], long double vectorVelocidad[3], long double vectorResultado[12])
{
    long double p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper;
    long double mu = 398600.4418e9;
    long double small = 1e-10;
    long double undefined = 999999.1;

    long double normaVectPosicion = Norma(vectorPosicion);
    long double normaVectVelocidad = Norma(vectorVelocidad);

    long double vHbar[3];
    cross(vectorPosicion, vectorVelocidad, vHbar);

    long double normaHbar = Norma(vHbar);

    long double saliadNewton[2];

    if ( normaHbar > small)
    {
        long double nbar[3];
        nbar[0] = -vHbar[1];
        nbar[1] =  vHbar[0];
        nbar[2] =   0.0;

        long double normaNbar = Norma(nbar);
        long double c1 = normaVectVelocidad*normaVectVelocidad - mu / normaVectPosicion;

        long double rdotv = dot(3, vectorPosicion, vectorVelocidad);
        long double ebar[3];
        for (int i=0; i<3; i++)
            ebar[i] = (c1*vectorPosicion[i] - rdotv*vectorVelocidad[i])/mu;

        ecc = Norma( ebar );

        long double sme = ( normaVectVelocidad*normaVectVelocidad*0.5  ) - ( mu /normaVectPosicion );

        if ( fabs( sme ) > small )
        {
            a = -mu  / (2.0 *sme);
        }
        else
        {
            a = 0.0;
        }

        p = normaHbar*normaHbar/mu;

        long double hk  = vHbar[2]/normaHbar;
        incl= acos( hk );

        char typeorbit[2]= "ei";
        if ( ecc < small )
        {
            //--------------  circular equatorial ---------------
            if  ((incl<small) || (fabs(incl-PI)<small))
            {
                typeorbit[0] = 'c';
                typeorbit[1] = 'e';
            }
            else
            {
                typeorbit[0] = 'c';
                typeorbit[1] = 'i';
            }

        }
        else
        {
            //elliptical, parabolic, hyperbolic equatorial --
            if  ((incl<small) || (fabs(incl-PI)<small))
            {
                typeorbit[0] = 'e';
                typeorbit[1] = 'e';
            }

        }

        //--------  find longitude of ascending node ------------
        long double temp = 0.0;
        if ( normaNbar > small )
        {
            temp = nbar[0] / normaNbar;
            if ( abs(temp) > 1.0  )
            {
                temp= Sing(temp);
            }
            omega= acos( temp );
            if ( nbar[1] < 0.0  )
            {
                omega= 2*PI - omega;
            }
        }
        else
        {
            omega= undefined;
        }

        //-------------- find argument of perigee ---------------

        if ( typeorbit[0] == 'e' && typeorbit[1] == 'i')
        {
            argp = angl( nbar,ebar);
            if ( ebar[2] < 0.0  )
                argp= 2*PI - argp;
        }
        else{
            argp= undefined;
        }



        //----------  find true anomaly at epoch    -------------
        if ( typeorbit[0] == 'e' )
        {
            nu =  angl( ebar, vectorPosicion);
            if ( rdotv < 0.0  )
                nu= 2*PI - nu;
        }
        else
            nu= undefined;

        //--  find argument of latitude - circular inclined -----
        if ( strcmp(typeorbit,"ci") == 0 )
        {
            arglat = angl( nbar, vectorPosicion );
            if ( vectorPosicion[2] < 0.0  )
                arglat= 2*PI - arglat;
            m = arglat;
        }
        else
            arglat= undefined;

        // -- find longitude of perigee - elliptical equatorial ----
        if  (( ecc>small ) && (( typeorbit[0] == 'e' && typeorbit[1] == 'e')) )
        {
            temp= ebar[0]/ecc;
            if ( fabs(temp) > 1.0  )
                temp= Sing(temp);
            lonper= acos( temp );
            if ( ebar[1] < 0.0  )
                lonper= 2*PI - lonper;
            if ( incl > PI/2 )
                lonper= 2*PI - lonper;
        }
        else
            lonper= undefined;

        //------ find true longitude - circular equatorial ------
        if  (( normaVectPosicion>small ) && ( typeorbit[0] == 'c' && typeorbit[1] == 'e') == 0 )
        {
            temp= vectorPosicion[0]/normaVectPosicion;
            if ( fabs(temp) > 1.0  )
                temp = Sing(temp);
            truelon = acos( temp );
            if ( vectorPosicion[1] < 0.0  )
                truelon= 2*PI - truelon;
            if ( incl > PI/2 )
                truelon= 2*PI - truelon;
            m = truelon;
        }
        else
            truelon= undefined;

        //---------- find mean anomaly for all orbits -----------
        if ( typeorbit[0] == 'e' ) {

            newtonnu(ecc,nu, saliadNewton );
            m = saliadNewton[1];
        }

    }
    else
    {
        p    = undefined;
        a    = undefined;
        ecc  = undefined;
        incl = undefined;
        omega= undefined;
        argp = undefined;
        nu   = undefined;
        m    = undefined;
        arglat = undefined;
        truelon= undefined;
        lonper = undefined;
    }

    vectorResultado[0] = p;
    vectorResultado[1] = a;
    vectorResultado[2] = ecc;
    vectorResultado[3] = incl;
    vectorResultado[4] = incl;
    vectorResultado[5] = omega;
    vectorResultado[6] = argp;
    vectorResultado[7] = nu;
    vectorResultado[8] = m;
    vectorResultado[9] = arglat;
    vectorResultado[10] = truelon;
    vectorResultado[11] = lonper;

    /*
    printf("p = %f \n", p);
    printf("a = %f \n", a);
    printf("ecc = %f \n", ecc);
    printf("incl = %f \n", incl);
    printf("omega = %f \n", omega);
    printf("argp = %f \n", argp);
    printf("nu = %f \n", nu);
    printf("m = %f \n", m);
    printf("arglat = %f \n", arglat);
    printf("truelon = %f \n", truelon);
    printf("lonper = %f \n", lonper);
*/
}

