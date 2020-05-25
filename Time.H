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

#ifndef TIME_H_
#define TIME_H_

#include <chrono>
#include <iostream>
#include <cmath>

static inline
long long int GetNow(void) {
  return std::chrono::duration_cast<std::chrono::microseconds>(
           std::chrono::system_clock::now().time_since_epoch()).count();
}

static constexpr double J2000 = 2451545.0;

static inline constexpr
double TimeToJD_minus_J2000(long long int time) {
  return (0.000001/86400) * time - (365*30+7+0.5);// + J2000;
}

  // Greenwich mean sidereal time
  // GMST in range [0..1):
static inline double JDMinusJ2000ToGMST(double jd_minus_j2000) {
    // Meeus
    // jd -= J2000;
    // const double T = jd  / 36525
    // theta0 = 280.46061837 + 360.98564736629*(JD-2451545.0) + 0.000387933*T*T - T*T*T/38710000.0
  const double x = ( ( (1.0/(38710000.0*36525.0*36525.0*36525.0*360.0)) * jd_minus_j2000 +
                       (0.000387933 / (36525.0*36525.0*360.0))) * jd_minus_j2000 +
                       (360.98564736629/360.0)) * jd_minus_j2000 +
                       (280.46061837/360.0);
  return x - floor(x);
}

  // GMST in range [0..1):
static inline
double TimeToGMST(long long int time) {
  return JDMinusJ2000ToGMST(TimeToJD_minus_J2000(time));
}

void DecomposeDay(int day,int &y,int &m,int &d);
    // day.. number of days since 1970.01.01

struct PrintTime {
  PrintTime(void) : time(GetNow()) {}
  PrintTime(long long int time) : time(time) {}
  const long long int time;
};

std::ostream &operator<<(std::ostream &o,const PrintTime &t);

#endif
