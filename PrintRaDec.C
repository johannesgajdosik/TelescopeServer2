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

std::ostream &operator<<(std::ostream &o,const PrintRaSeconds &x) {
  char buff[9];
  char *p = buff + sizeof(buff);
  unsigned int i = x.ra_secs / 10;
  *--p = '\0';
  *--p = '0' + x.ra_secs - 10*i;
  unsigned int j = i / 6;
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
};

std::ostream &operator<<(std::ostream &o,const PrintDecSeconds &x) {
  char buff[10];
  char *p = buff + sizeof(buff);
  unsigned int i = x.dec_secs / 10;
  *--p = '\0';
  *--p = '0' + x.dec_secs - 10*i;
  unsigned int j = i / 6;
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
}

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
}

