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
#include "AcceptPort.H"
#include "PrintRaDec.H"
#include "Time.H"

#include <boost/asio.hpp>
#include <boost/bind.hpp>

#include <thread>
#include <chrono>
#include <iostream>


static Telescope::Ptr telescope;

static void HandleMsgGoto(unsigned int ra_int,int dec_int) {
//  std::cout << "HandleMsgGoto(" << PrintRaInt(ra_int) << ',' << PrintDecInt(dec_int) << ')'
//            << std::endl;
  if (telescope) telescope->gotoPosition(ra_int,dec_int);
}

static void HandleMsgGuide(int d_ra_micros,int d_dec_micros) {
//  std::cout << "HandleMsgGuide(" << d_ra_micros << ',' << d_dec_micros << ')'
//            << std::endl;
  if (telescope) telescope->guide(d_ra_micros,d_dec_micros);
}

static void HandleMsgMove(short int horz,short int vert,unsigned int micros) {
//  std::cout << "HandleMsgMove(" << horz << ',' << vert << ';' << micros << ')'
//            << std::endl;
  if (telescope) telescope->move(horz,vert,micros);
}

static int HandleStellariumTelescopeProtocolServerMsg(
             const void *data,int size,TcpIpConnection &c) {
  if (size < 4) return 0;
  const int length = reinterpret_cast<const unsigned short*>(data)[0];
  const unsigned short type = reinterpret_cast<const unsigned short*>(data)[1];
  if (size < length) return 0;
  switch (type) {
    case 0:
      if (length < 20) return -1;
      HandleMsgGoto(
        reinterpret_cast<const unsigned int*>(data)[3],
        reinterpret_cast<const int*>(data)[4]);
      break;
    case 1:
      if (length < 20) return -1;
      HandleMsgGuide(
        reinterpret_cast<const int*>(data)[3],
        reinterpret_cast<const int*>(data)[4]);
      break;
    case 2:
      if (length < 20) return -1;
      HandleMsgMove(
        reinterpret_cast<const short*>(data)[6],
        reinterpret_cast<const short*>(data)[7],
        reinterpret_cast<const unsigned int*>(data)[4]);
    default:
        // just ignore
      break;
  }
  return length;
}

static
void BroadcastPosition(unsigned int ra_int,int dec_int,AcceptPort &server) {
  unsigned char msg[24];
  *reinterpret_cast<unsigned short*>(&msg) = 24;
  *reinterpret_cast<unsigned short*>(&msg[2]) = 0;
  *reinterpret_cast<unsigned long long*>(&msg[4]) = 0;
  *reinterpret_cast<unsigned int*>(&msg[12]) = ra_int;
  *reinterpret_cast<int*>(&msg[16]) = dec_int;
  *reinterpret_cast<int*>(&msg[20]) = 0;
  server.broadcast(msg,24);
}


static bool continue_looping = true;

static
void SignalHandler(boost::asio::signal_set &signals,
                   const boost::system::error_code &error,
                   int signal_number) {
  if (error) {
      // handler was cancelled?
  } else {
    signals.async_wait(boost::bind(SignalHandler,boost::ref(signals),_1,_2));
    std::cout << PrintTime() << " "
                 "signal " << signal_number << " received" << std::endl;
    continue_looping = false;
    throw nullptr;
  }
}

int main(int argc,char *argv[]) {
  std::signal(SIGPIPE,SIG_IGN);
  {
    boost::asio::io_context io_context;

    boost::asio::signal_set signals(io_context,SIGINT,SIGTERM);
    signals.async_wait(boost::bind(SignalHandler,boost::ref(signals),_1,_2));

    int port;
    if (argc < 3 || !(std::istringstream(argv[1])>>port) ||
        port <= 0 || port >= 0x10000) {
      std::cerr << "Usage: " << argv[0] << " <port> <telescope_type:args>" << std::endl;
      return 1;
    }

    AcceptPort server(io_context,port,
                      &HandleStellariumTelescopeProtocolServerMsg);

      telescope = Telescope::Create(argv[2],
                                    [&server](bool opened_closed) {
                                      if (opened_closed) server.start();
                                      else server.stop();
                                    },
                                    [&server](unsigned int ra_int,int dec_int) {
                                      BroadcastPosition(ra_int,dec_int,server);
                                    },
                                    io_context);
      if (telescope) {
//        try {
          io_context.run();
//        } catch (std::exception &e) {
//          std::cerr << PrintTime() << " "
//                       "Exception: " << e.what() << std::endl;
//          abort();
//        } catch (std::nullptr_t &) {
//          if (!continue_looping) break;
//          std::cerr << PrintTime() << " "
//                       "unrecoverable error" << std::endl;
//        }
        telescope.reset();
      }
      // TODO: the serial connection is destroyed after the telescope.
      // methods of the telescope object could be called after destruction.
      // io_context is destroyed here (end of scope), this will also release the serial connection.
      // TODO: further investigation
  }
  std::cout << PrintTime() << " "
               "main: bye." << std::endl;
  return 0;
}

