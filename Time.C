/*
 * Author and Copyright
 * Johannes Gajdosik, 2020
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include "Time.H"

void DecomposeDay(int day,int &y,int &m,int &d) {
    // day.. number of days since 1970.01.01
  day += (5*(400*365+97) -(30*365+7) -(31+29));
    // number of days since 0000.03.01
  y = day / (400*365+97);
  day -= (400*365+97) * y;
  if (day < 0) {day += (400*365+97);--y;}
  y *= 400;
  {
    const int c = day / (100*365+24);
    y += 100 * c;
    if (c > 3) {
      m = 1;
      d = 28;
      return;
    }
    day -= (100*365+24) * c;
  }
  {
    const int h = day / (4*365+1);
    day -= (4*365+1) * h;
    y += 4 * h;
  }
  {
    const int h = day / 365;
    y += h;
    if (h > 3) {
      m = 1;
      d = 28;
      return;
    }
    day -= 365 * h;
  }
  if (day < 31+30+31+30+31+31)
    if (day < 31+30+31)
      if (day < 31)
        {m = 2; d = day;} // mar
      else
       if (day < 31+30)
         {m = 3; d = day-31;} // apr
       else
         {m = 4; d = day-(31+30);} // may
    else
      if (day < 31+30+31+30)
        {m = 5; d = day-(31+30+31);} // jun
      else
        if (day < 31+30+31+30+31)
          {m = 6; d = day-(31+30+31+30);} // jul
        else
          {m = 7; d = day-(31+30+31+30+31);} // aug
  else
    if (day < 31+30+31+30+31+31+30+31+30)
      if (day < 31+30+31+30+31+31+30)
        {m = 8; d = day-(31+30+31+30+31+31);} // sep
      else
        if (day < 31+30+31+30+31+31+30+31)
          {m = 9; d = day-(31+30+31+30+31+31+30);} // oct
        else
          {m = 10; d = day-(31+30+31+30+31+31+30+31);} // nov
    else
      if (day < 31+30+31+30+31+31+30+31+30+31)
        {m = 11; d = day-(31+30+31+30+31+31+30+31+30);} // dec
      else {
        if (day < 31+30+31+30+31+31+30+31+30+31+31)
          {m = 0; d = day-(31+30+31+30+31+31+30+31+30+31);} // jan
        else
          {m = 1; d = day-(31+30+31+30+31+31+30+31+30+31+31);} // feb
        ++y;
      }
}

std::ostream &operator<<(std::ostream &o,const PrintTime &t) {
  char buff[24];
  char *p = buff+sizeof(buff);
  *--p = '\0';
  const int day = (t.time / 86400000000LL) - ((t.time < 0) ? 1 : 0);
  {
    long long int x = t.time - (day*86400000000LL);
    x /= 1000;
    {const long long int y = x / 10; *--p = '0' + (char)(x-y*10); x = y;}
    {const long long int y = x / 10; *--p = '0' + (char)(x-y*10); x = y;}
    {const long long int y = x / 10; *--p = '0' + (char)(x-y*10); x = y;}
    *--p = '.';
    {const long long int y = x / 10; *--p = '0' + (char)(x-y*10); x = y;}
    {const long long int y = x /  6; *--p = '0' + (char)(x-y* 6); x = y;}
    *--p = ':';
    {const long long int y = x / 10; *--p = '0' + (char)(x-y*10); x = y;}
    {const long long int y = x /  6; *--p = '0' + (char)(x-y* 6); x = y;}
    *--p = ':';
    {const long long int y = x / 10; *--p = '0' + (char)(x-y*10);
     *--p = '0' + (char)y;}
  }
  *--p = '-';
  int y,m,d;
  DecomposeDay(day,y,m,d);d++;m++;
  {const long long int x = d / 10; *--p = '0' + (char)(d-x*10);
   *--p = '0' + (char)x;}
  *--p = '.';
  {const long long int x = m / 10; *--p = '0' + (char)(m-x*10);
   *--p = '0' + (char)x;}
  *--p = '.';
  {const long long int x = y / 10; *--p = '0' + (char)(y-x*10); y = x;}
  {const long long int x = y / 10; *--p = '0' + (char)(y-x*10); y = x;}
  {const long long int x = y / 10; *--p = '0' + (char)(y-x*10);
   *--p = '0' + (char)x;}
  return o << p;
}

