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

#ifndef TELESCOPE_H_
#define TELESCOPE_H_

#include "IntrusiveBase.H"

#include <boost/intrusive_ptr.hpp>

namespace boost {namespace asio {class io_context;}}

#include <functional>
#include <string>

class Telescope : public IntrusiveBase {
public:
  typedef boost::intrusive_ptr<Telescope> Ptr;
  typedef std::function<void(unsigned int ra_int,int dec_int)
                       > PositionFunction;
  typedef std::function<void(bool opened_closed)
                       > OpenedClosedFunction;
  typedef std::function<Ptr(const std::string &args,
                            OpenedClosedFunction &&opened_closed,
                            PositionFunction &&announce_position,
                            boost::asio::io_context &io_context)
                       > CreationFunction;
  static Ptr Create(const std::string &type_and_args,
                    OpenedClosedFunction &&opened_closed,
                    PositionFunction &&announce_position,
                    boost::asio::io_context &io_context);
  static bool RegisterCreationFunction(const std::string &type,
                                       CreationFunction &&f);
  virtual ~Telescope(void) {}
  virtual void gotoPosition(unsigned int ra_int,int dec_int) = 0;
  virtual void guide(int d_ra_micros,int d_dec_micros) = 0;
  virtual void move(short int horz,short int vert,
                    unsigned int validity_micros) = 0;
private:
  virtual bool initializationOk(void) const = 0;
protected:
  Telescope(OpenedClosedFunction &&opened_closed,
            PositionFunction &&announce_position)
    :opened_closed(std::move(opened_closed)),
     announce_position(std::move(announce_position)) {}
  const OpenedClosedFunction opened_closed;
  const PositionFunction announce_position;
};

#endif
