
/**
    Diferencia de tiempos

    Last modified:   2015/08/12   M. Mahooti
*/

void timeDiff( double UT1_UTC,  double TAI_UTC,  double diferenciaTiempos[5])
{
     double TT_TAI  = +32.184;
     double GPS_TAI = -19.0;

     double UT1_TAI = UT1_UTC-TAI_UTC;
     double UTC_TAI = -TAI_UTC;
     double UTC_GPS = UTC_TAI-GPS_TAI;
     double UT1_GPS = UT1_TAI-GPS_TAI;
     double TT_UTC  = TT_TAI-UTC_TAI;
     double GPS_UTC = GPS_TAI-UTC_TAI;

    diferenciaTiempos[0] = UT1_TAI;
    diferenciaTiempos[1] = UTC_GPS;
    diferenciaTiempos[2] = UT1_GPS;
    diferenciaTiempos[3] = TT_UTC;
    diferenciaTiempos[4] = GPS_UTC;

}
