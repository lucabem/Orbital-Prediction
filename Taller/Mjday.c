
#include "Mjday.h"

long double Mjday(int year, int month,int day, int hour, int min, long double sec)
{

    int y = year;
    int m = month;
    long double b = 0.0;
    long double c = 0.0;
    long double a = 0.0;
    if (m <=2)
    {
        y = y-1;
        m = m+12;
    }
    if(y<0)
    {
        c = -0.75;
    }

    // check for valid calendar date
    if (year < 1582)
    {

    }
    else if (year > 1582)
    {
        a = fix(y/100);
        b = 2 - a + floor(a/4);

    }
    else if ( month < 10)
    {

    }
    else if (month > 10)
    {
        a = fix(y/100);
        b = 2 - a + floor(a/4);
    }
    else if (day <= 4)
    {

    }
    else if (day > 14)
    {
        a = fix(y/100);
        b = 2 - a + floor(a/4);
    }
    else
    {
        return 0.0;
    }


    long double jd = fix(365.25*y+c) + fix(30.6001 * (m + 1));

    jd = jd + day + b + 1720994.5;

    jd = jd + (hour*1.0+min*1.0/60+sec/3600)/24;

    long double Mjd = jd - 2400000.5;

    return Mjd;
}
