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

#include "PrintRaDec.H"

std::ostream &operator<<(std::ostream &o,const PrintRaMilliseconds &x) {
  char buff[13];
  char *p = buff + sizeof(buff);
  unsigned int j = x.ra_millis / 10;
  *--p = '\0';
  *--p = '0' + x.ra_millis - 10*j;
  unsigned int i = j / 10;
  *--p = '0' + j - 10*i;
  j = i/10;
  *--p = '0' + i - 10*j;
  *--p = 's';
  i = j/10;
  *--p = '0' + j - 10*i;
  j = i/6;
  *--p = '0' + i - 6*j;
  *--p = 'm';
  i = j/10;
  *--p = '0' + j - 10*i;
  j = i/6;
  *--p = '0' + i - 6*j;
  *--p = 'h';
  i = j/10;
  *--p = '0' + j - 10*i;
  *--p = '0' + i;
  if (i>2) abort();
  return o << buff;
/*
  unsigned int hours = x.ra_millis/1000;
  const unsigned int milliseconds = x.ra_millis - hours*1000;
  unsigned int i = hours;
  hours = i/60;
  const unsigned int seconds = i - hours*60;
  i = hours;
  hours = i/60;
  const unsigned int minutes = i - hours*60;
  o << std::setw(2) << std::setfill('0') << hours << 'h'
    << std::setw(2) << std::setfill('0') << minutes << 'm'
    << std::setw(2) << std::setfill('0') << seconds << 's'
    << std::setw(3) << std::setfill('0') << milliseconds;
  return o;
*/
};

std::ostream &operator<<(std::ostream &o,const PrintDecCentiseconds &x) {
  char buff[13];
  char *p = buff + sizeof(buff);
  unsigned int i = x.dec_centis / 10;
  *--p = '\0';
  *--p = '0' + x.dec_centis - 10*i;
  unsigned int j = i / 10;
  *--p = '0' + i - 10*j;
  *--p = 's';
  i = j/10;
  *--p = '0' + j - 10*i;
  j = i/6;
  *--p = '0' + i - 6*j;
  *--p = 'm';
  i = j/10;
  *--p = '0' + j - 10*i;
  j = i/6;
  *--p = '0' + i - 6*j;
  *--p = 'd';
  i = j/10;
  *--p = '0' + j - 10*i;
  *--p = '0' + i;
  if (i>9) abort();
  buff[0] = (x.dec_sign?'-':'+');
  return o << buff;
/*

  unsigned int degrees = x.dec_centis/100;
  const unsigned int centiseconds = x.dec_centis - degrees*100;
  unsigned int i = degrees;
  degrees = i/60;
  const unsigned int seconds = i - degrees*60;
  i = degrees;
  degrees = i/60;
  const unsigned int minutes = i - degrees*60;
  o << (x.dec_sign?'-':'+')
    << std::setw(2) << std::setfill('0') << degrees << 'd'
    << std::setw(2) << std::setfill('0') << minutes << 'm'
    << std::setw(2) << std::setfill('0') << seconds << 's'
    << std::setw(2) << std::setfill('0') << centiseconds;
  return o;
*/
}

