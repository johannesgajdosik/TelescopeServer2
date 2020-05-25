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

#include "Connection.H"
#include "PrintRaDec.H"
#include "Time.H"

#include <boost/asio.hpp>
#include <boost/bind.hpp>

#include <thread>
#include <chrono>
#include <iostream>

#include <fcntl.h>
#include <linux/joystick.h>

static TcpIpConnection::Ptr client;

typedef Connection<boost::asio::posix::stream_descriptor> PosixConnection;
static PosixConnection::Ptr joystick_client;

static void HandleMsgCurrentPosition(unsigned int ra_int,int dec_int) {
//  std::cout << PrintTime() << " "
//               "HandleMsgCurrentPosition(" << PrintRaInt(ra_int) << ','
//            << PrintDecInt(dec_int) << ')' << std::endl;
}

static void SendMsgGoto(unsigned int ra_int,int dec_int) {
  if (!client) return;
  unsigned char msg[24];
  *reinterpret_cast<unsigned short*>(&msg) = 24;
  *reinterpret_cast<unsigned short*>(&msg[2]) = 0;
  *reinterpret_cast<unsigned long long*>(&msg[4]) = 0;
  *reinterpret_cast<unsigned int*>(&msg[12]) = ra_int;
  *reinterpret_cast<int*>(&msg[16]) = dec_int;
  *reinterpret_cast<int*>(&msg[20]) = 0;
  client->write(msg,24);
}

static void SendMsgMove(short int horz,short int vert,unsigned int micros) {
  std::cout << PrintTime() << " "
               "SendMsgMove(" << horz << ',' << vert << ';' << micros << ')' << std::endl;
  if (!client) return;
  unsigned char msg[20];
  *reinterpret_cast<unsigned short*>(&msg) = 20;
  *reinterpret_cast<unsigned short*>(&msg[2]) = 2;
  *reinterpret_cast<unsigned long long*>(&msg[4]) = 0;
  *reinterpret_cast<short*>(&msg[12]) = horz;
  *reinterpret_cast<short*>(&msg[14]) = vert;
  *reinterpret_cast<unsigned int*>(&msg[16]) = micros;
  client->write(msg,20);
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
      HandleMsgCurrentPosition(
        reinterpret_cast<const unsigned int*>(data)[3],
        reinterpret_cast<const int*>(data)[4]);
      break;
    default:
        // just ignore
      break;
  }
  return length;
}



static std::unique_ptr<boost::asio::deadline_timer> joystick_deadline;

static short int value_x = 0;
static short int value_y = 0;
static unsigned int buttons = 0;


static void SendMsgMove(void) {
  SendMsgMove(value_x,value_y,1000000);
  if (value_x == 0 && value_y == 0) return;
  if (!joystick_deadline) return;
  joystick_deadline->expires_from_now(boost::posix_time::microseconds(500000));
  joystick_deadline->async_wait(
                       [](const boost::system::error_code &e) {
                         if (e) return;
                         SendMsgMove();
                       });
}

static void HandleJoystickEvent(const struct js_event &event) {
  switch (event.type) {
    case JS_EVENT_BUTTON:
      if ((unsigned int)event.number < 32) {
        if (event.value != 0) {
          buttons |= 1u << ((unsigned int)event.number);
        } else {
          buttons &= ~(1u << ((unsigned int)event.number));
        }
      }
      std::cout << PrintTime() << " "
                   "button " << (unsigned int)event.number << ": " << event.value << std::endl;
      break;
    case JS_EVENT_AXIS:
      std::cout << PrintTime() << " "
                   "Axis(" << (unsigned int)event.number << "): " << event.value << std::endl;
      switch ((unsigned int)event.number) {
        case 0:
        case 3:
        case 4:
            // right is +
          std::cout << PrintTime() << " "
                       "x: " << event.value << std::endl;
          if (event.value == 0) {
            value_x = 0;
            SendMsgMove();
          } else {
            int x = buttons & 15;
            if (x) {
              if (x > 9) x = 9;
              value_x = (x*event.value)/9;
              SendMsgMove();
            }
          }
          SendMsgMove();
          break;
        case 1:
        case 2:
        case 5:
            // down is +
          std::cout << PrintTime() << " "
                       "y: " << event.value << std::endl;
          if (event.value == 0) {
            value_y = 0;
            SendMsgMove();
          } else {
            short int x = buttons & 15;
            if (x) {
              if (x > 9) x = 9;
              value_y = (x*event.value)/9;
              SendMsgMove();
            }
          }
          break;
      }
      break;
  }
}

static int HandleJoystickMessage(const void *data,int size,PosixConnection&) {
  int rval = 0;
  while (size >= (int)sizeof(struct js_event)) {
    HandleJoystickEvent(*((const struct js_event*)data));
    data = ((struct js_event*)data) + 1;
    size -= sizeof(struct js_event);
    rval += sizeof(struct js_event);
  }
  return rval;
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

  int port;
  if (argc < 4 || !(std::istringstream(argv[2])>>port) ||
      port <= 0 || port >= 0x10000) {
    std::cerr << "Usage: " << argv[0] << " <host> <port> </dev/input/js0 or similar>\n\n"
              << "Hint: The 4 buttons define a bitfield with values 0..15 and control the movement speed. "
                 "Values >9 are cropped to 9. Only while at least one button is pressed the scope will move at all. "
                 "This is for avoiding unvoluntary moves during exposure. "
                 "With analog axes you can additionally fine tune the movement speed."
              << std::endl;
    return 1;
  }

  int joy_fd = open(argv[3],O_RDONLY);
  if (joy_fd < 0) {
    std::cerr << "open(" << argv[3] << ") failed: " << strerror(errno) << std::endl;
    return 1;
  }
  if (fcntl(joy_fd,F_SETFL,O_NONBLOCK) < 0) {
    std::cerr << "fcntl(O_NONBLOCK) failed: " << strerror(errno) << std::endl;
    close(joy_fd);
    return 1;
  }

  {
    boost::asio::io_context io_context;
    joystick_deadline = std::make_unique<boost::asio::deadline_timer>(io_context);

    boost::asio::signal_set signals(io_context,SIGINT,SIGTERM);
    signals.async_wait(boost::bind(SignalHandler,boost::ref(signals),_1,_2));

    boost::asio::posix::stream_descriptor joystick(io_context,joy_fd);

    joystick_client
      = PosixConnection::Create(std::move(joystick),
                                HandleJoystickMessage,
                                PosixConnection::ErrorHandler());

    while (continue_looping) {
      boost::asio::ip::tcp::resolver resolver(io_context);
      boost::asio::ip::tcp::socket socket(io_context);

      resolver.async_resolve(
        boost::asio::ip::tcp::resolver::query(argv[1],argv[2]),
        [s=&socket]
            (const boost::system::error_code &ec,
             boost::asio::ip::tcp::resolver::results_type endpoints) {
          if (ec) {
            std::cerr << PrintTime() << " "
                         "resolving failed: " << ec.message() << std::endl;
            throw nullptr;
          }
          std::cerr << PrintTime() << " "
                       "start connecting..." << std::endl;
          boost::asio::async_connect(
            *s,endpoints,
            [s](const boost::system::error_code &ec,
                     const boost::asio::ip::tcp::endpoint &endpoint) {
              if (ec) {
                std::cerr << PrintTime() << " "
                             "connect failed: " << ec.message() << std::endl;
                throw nullptr;
              }
                // disable Nagler algorithm
              s->set_option(boost::asio::ip::tcp::no_delay(true));

              std::cerr << PrintTime() << " "
                           "connect ok" << std::endl;
              client
                = TcpIpConnection::Create(
                    std::move(*s),
                    HandleStellariumTelescopeProtocolServerMsg,
                    [](const boost::system::error_code &error,
                       TcpIpConnection &c) {
                      throw nullptr;
                    });
            });
        });
      try {
        io_context.run();
      } catch (std::exception &e) {
        std::cerr << PrintTime() << " "
                     "Exception: " << e.what() << "\n";
      } catch (std::nullptr_t &) {
        if (!continue_looping) break;
        std::cerr << PrintTime() << " "
                     "unrecoverable error" << std::endl;
      }
      client.reset();
      std::this_thread::sleep_for(std::chrono::microseconds(1000000));
      joystick_deadline.reset();
      client.reset();
      joystick_client.reset();
    }
  }
  close(joy_fd);
  std::cout << PrintTime() << " "
               "main: bye." << std::endl;
  return 0;
}
