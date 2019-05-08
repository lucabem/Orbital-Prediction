#include "NewtonNu.h"
#include <math.h>
#include "MatlabUtilities.h"
#include "Constantes.h"

void newtonnu( double ecc,  double nu,  double salida[2])
{
    double e0= 999999.9;
    double m = 999999.9;
    double small = 0.00000001;

// --------------------------- circular ------------------------
    if ( fabs( ecc ) < small  )
    {
        m = nu;
        e0= nu;
    }
    else
    {
        // ---------------------- elliptical -----------------------
        if ( ecc < 1.0-small  )
        {
            double sine= ( sqrt( 1.0 -ecc*ecc ) * sin(nu) ) / ( 1.0 +ecc*cos(nu) );
            double cose= ( ecc + cos(nu) ) / ( 1.0  + ecc*cos(nu) );
            e0  = atan2( sine,cose );
            m   = e0 - ecc*sin(e0);
        }
        else
        {
            // -------------------- hyperbolic  --------------------
            if ( ecc > 1.0 + small  )
            {
                if ((ecc > 1.0 ) && (fabs(nu)+0.00001 < PI-acos(1.0 /ecc)))
                {
                    double sine= ( sqrt( ecc*ecc-1.0  ) * sin(nu) ) / ( 1.0  + ecc*cos(nu) );
                    e0  = asinh( sine );
                    m   = ecc*sinh(e0) - e0;
                }
            }
            else
            {
                // ----------------- parabolic ---------------------
                if ( fabs(nu) < 168.0*PI/180.0  )
                {
                    e0= tan( nu*0.5  );
                    m = e0 + (e0*e0*e0)/3.0;
                }
            }
        }

    }

    if ( ecc < 1.0  )
    {
        m = fmod( m,2.0 *PI );
        if ( m < 0.0  )
            m= m + 2.0 *PI;
        e0 = fmod( e0,2.0 *PI );
    }

    salida[0] = e0;
    salida[1] = m;
}
