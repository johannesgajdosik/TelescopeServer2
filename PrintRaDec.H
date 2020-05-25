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

#ifndef PRINT_RA_DEC_H
#define PRINT_RA_DEC_H

#include <iostream>
#include <cmath> // floor

struct PrintRaMilliseconds {
  const unsigned int ra_millis;
  PrintRaMilliseconds(unsigned int ra_millis) : ra_millis(ra_millis) {}
};

struct PrintRaInt : public PrintRaMilliseconds {
  PrintRaInt(unsigned int ra_int) : PrintRaMilliseconds( (ra_int*86400000ULL + 0x80000000ULL)>>32 ) {}
};

std::ostream &operator<<(std::ostream &o,const PrintRaMilliseconds &x);

static inline
unsigned int RadiantToRaMillis(double ra_rad) {
  ra_rad *= (1.0/(2*M_PI));
    // ((ra_rad-floor(ra_rad)) * 86400000.0 + 0.5) might actually be >= 86400000.
    // In this case conversion to (unsigned int) would be undefined,
    // so I explicitely convert to unsigned long long int,
    // which is converted to unsigned int modulo (1<<32)
  return (unsigned long long int)((ra_rad-floor(ra_rad)) * 86400000.0 + 0.5);
}

static inline
unsigned int RadiantToUnsignedInt(double ra_rad) {
  ra_rad *= (1.0/(2*M_PI));
    // ((ra_rad-floor(ra_rad)) * 86400000.0 + 0.5) might actually be >= 86400000.
    // In this case conversion to (unsigned int) would be undefined,
    // so I explicitely convert to unsigned long long int,
    // which is converted to unsigned int modulo (1<<32)
  return (unsigned long long int)((ra_rad-floor(ra_rad)) * 4294967296.0 + 0.5);
}

struct PrintRaRad : public PrintRaMilliseconds {
  PrintRaRad(double ra_rad) : PrintRaMilliseconds(RadiantToRaMillis(ra_rad)) {}
};

static inline
unsigned int RadiantToDecCentis(double dec_rad) {
  dec_rad *= (1.0/(2*M_PI));
  return (unsigned long long int)((dec_rad-floor(dec_rad)) * (360*3600*100.0) + 0.5);
}

struct PrintDecCentiseconds {
  const bool dec_sign;
  const unsigned int dec_centis;
  PrintDecCentiseconds(bool dec_sign,unsigned int dec_centis)
      : dec_sign(dec_sign),dec_centis(dec_centis) {
    if (dec_centis > 90*3600*100) abort();
  }
protected:
  static PrintDecCentiseconds Create(double dec_rad) {
    const bool dec_sign = (dec_rad < 0.0);
    return PrintDecCentiseconds(dec_sign,RadiantToDecCentis(dec_sign ? (-dec_rad) : dec_rad));
  }
};

std::ostream &operator<<(std::ostream &o,const PrintDecCentiseconds &x);

struct PrintDecInt : public PrintDecCentiseconds {
  PrintDecInt(int dec_int)
    : PrintDecCentiseconds((dec_int<0),
                            // -dec_int might cause int overflow, -(unsigned int)dec_int is ok.
                           (unsigned int)( ( ( ((dec_int<0) ? (-(unsigned int)dec_int) : (unsigned int)dec_int)
                             * (360*3600*100ULL) ) + 0x80000000ULL ) >> 32 ) ) {}
};

struct PrintDecRad : public PrintDecCentiseconds {
  PrintDecRad(double dec_rad) : PrintDecCentiseconds(PrintDecCentiseconds::Create(dec_rad)) {}
};


#endif
