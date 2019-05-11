

void timeDiff( long double UT1_UTC,  long double TAI_UTC,  long double diferenciaTiempos[5])
{
     long double TT_TAI  = +32.184;
     long double GPS_TAI = -19.0;
//     long double TT_GPS  =  TT_TAI-GPS_TAI;
//     long double TAI_GPS = -GPS_TAI;
     long double UT1_TAI = UT1_UTC-TAI_UTC;
     long double UTC_TAI = -TAI_UTC;
     long double UTC_GPS = UTC_TAI-GPS_TAI;
     long double UT1_GPS = UT1_TAI-GPS_TAI;
     long double TT_UTC  = TT_TAI-UTC_TAI;
     long double GPS_UTC = GPS_TAI-UTC_TAI;

    diferenciaTiempos[0] = UT1_TAI;
    diferenciaTiempos[1] = UTC_GPS;
    diferenciaTiempos[2] = UT1_GPS;
    diferenciaTiempos[3] = TT_UTC;
    diferenciaTiempos[4] = GPS_UTC;

}
