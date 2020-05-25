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

#include "Telescope.H"

#include <boost/asio.hpp>

#include <iostream>

class TelescopeDummy : public Telescope {
public:
  static Ptr Create(const std::string &args,
                    OpenedClosedFunction &&opened_closed,
                    PositionFunction &&announce_position,
                    boost::asio::io_context &io_context) {
    return new TelescopeDummy(args,
                              std::move(opened_closed),
                              std::move(announce_position),
                              io_context);
  }
private:
  static int Speed(int x) {
      // TODO: something more sophisticated, like handling maximimum speed
    if (x >= 0) return (x+1)/2;
    return (x-1)/2;
  }
  void handlePeriodic(void) {
    pos_ra += Speed(goto_pos_ra - pos_ra);
    pos_dec += Speed(goto_pos_dec - pos_dec);
    announce_position(pos_ra,pos_dec);
  }
  void setupTimer(void) {
    timer.async_wait(
      [this](const boost::system::error_code &e) {
        if (e) return;
        handlePeriodic();
        timer.expires_at(timer.expires_at()
                           + boost::posix_time::microseconds(250000));
        setupTimer();
      });
  }
  bool initializationOk(void) const override {return true;}
  TelescopeDummy(const std::string &args,
                 OpenedClosedFunction &&opened_closed,
                 PositionFunction &&announce_position,
                 boost::asio::io_context &io_context)
      : Telescope(std::move(opened_closed),std::move(announce_position)),
        timer(io_context) {
    timer.expires_from_now(boost::posix_time::microseconds(0));
    setupTimer();
    opened_closed(true);
  }
  ~TelescopeDummy(void) {
    opened_closed(false);
  }
  void gotoPosition(unsigned int ra_int,int dec_int) override {
    std::cout << "TelescopeDummy::gotoPosition("
              << ra_int << ',' << dec_int << ')' << std::endl;
    goto_pos_ra = ra_int;
    goto_pos_dec = dec_int;
  }
  void guide(int d_ra_micros,int d_dec_micros) override {
  }
  void move(short int horz,short int vert,unsigned int micros) override {
  }
  boost::asio::deadline_timer timer;
  unsigned int pos_ra = 0;
  int pos_dec = 0;
  unsigned int goto_pos_ra = 0;
  int goto_pos_dec = 0;
};

static const bool telescope_dummy_registered
  = Telescope::RegisterCreationFunction("Dummy",TelescopeDummy::Create);

