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
#include "Telescope.H"

#include <boost/asio.hpp>

#include <iostream>
#include <string>

static TcpIpConnection::Ptr client;

static int ParseCoordinateSystem(const char str[]) {
  if (!str) return -1;
  if (0 == strcmp(str,"J2000")) return 0;
  if (0 == strcmp(str,"equatorial")) return 1;
  if (0 == strcmp(str,"azimut")) return 2;
  if (0 == strcmp(str,"land")) return 3;
  return -1;
}

const char *CoordinateSystemToString(Telescope::CoordinateSystem coordinate_system) {
  switch (coordinate_system) {
    case Telescope::CoordinateSystem::J2000: return "J2000";
    case Telescope::CoordinateSystem::equatorial: return "equatorial";
    case Telescope::CoordinateSystem::azimut: return "azimut";
    case Telescope::CoordinateSystem::land: return "land";
  }
  return "undefined";
}

static void HandleMsgCurrentPosition(unsigned int ra_int,int dec_int,
                                     Telescope::CoordinateSystem coordinate_system) {
  std::cout << "HandleMsgCurrentPosition(" << ra_int << ',' << dec_int << ','
            << CoordinateSystemToString(coordinate_system) << ')'
            << std::endl;
}

static void SendMsgGoto(unsigned int ra_int,int dec_int,
                        Telescope::CoordinateSystem coordinate_system) {
  static constexpr int msg_length = 28;
  unsigned char msg[msg_length];
  *reinterpret_cast<unsigned short*>(&msg) = msg_length;
  *reinterpret_cast<unsigned short*>(&msg[2]) = 0;
  *reinterpret_cast<unsigned long long*>(&msg[4]) = 0; // time
  *reinterpret_cast<unsigned int*>(&msg[12]) = ra_int;
  *reinterpret_cast<int*>(&msg[16]) = dec_int;
  *reinterpret_cast<int*>(&msg[20]) = 0; // status
  *reinterpret_cast<int*>(&msg[24]) = static_cast<int>(coordinate_system);
  client->write(msg,msg_length);
}

static int HandleStellariumTelescopeProtocolServerMsg(
             const void *data,int size,TcpIpConnection &c) {
  if (size < 4) return 0;
  const int msg_length = reinterpret_cast<const unsigned short*>(data)[0];
  const unsigned short type = reinterpret_cast<const unsigned short*>(data)[1];
  if (size < msg_length) return 0;
  switch (type) {
    case 0:
      if (msg_length < 20) return -1;
      {
        const unsigned int ra_int = reinterpret_cast<const unsigned int*>(data)[3];
        const int dec_int = reinterpret_cast<const int*>(data)[4];
        Telescope::CoordinateSystem coordinate_system
          = Telescope::CoordinateSystem::J2000;
        if (msg_length >= 28) {
          coordinate_system
            = reinterpret_cast<const Telescope::CoordinateSystem*>(data)[6];
        }
        HandleMsgCurrentPosition(ra_int,dec_int,coordinate_system);
      }
      break;
    default:
        // just ignore
      break;
  }
  return msg_length;
}

int ParseUInt(char *&str) {
  if (!str) return -1;
  if (!*str) return -2;
  if (*str < '0' || *str > '9') return -3;
  int rval = *str++ - '0';
  while (*str >= '0' && *str <= '9') {
    rval = 10*rval + (*str++ - '0');
  }
  return rval;
}

static
int ScanDorHMS(char *str,char d_or_h) {
  if (!str) {
std::cout << "E1" << std::endl;
    return std::numeric_limits<int>::min();
  }
  int sign = 1;
  if (*str == '+') str++;
  else if (*str == '-') {sign = -1;str++;}
  if (d_or_h=='h' && sign < 0) {
std::cout << "E2" << std::endl;
    return std::numeric_limits<int>::min();
  }
  const int h = ParseUInt(str);
  if (h < 0) return std::numeric_limits<int>::min();
  if ((d_or_h=='h' && h>=24) || (d_or_h=='d' && h>=90)) {
std::cout << "E3" << std::endl;
    return std::numeric_limits<int>::min();
  }
  if (*str == '\0') return sign*h*3600;
  if (*str != d_or_h) {
std::cout << "E4" << std::endl;
    return std::numeric_limits<int>::min();
  }
  str++;
  if (*str == '\0') return sign*h*3600;
  const int m = ParseUInt(str);
  if (m < 0) {
std::cout << "E5" << std::endl;
    return std::numeric_limits<int>::min();
  }
  if (h==90 && m > 0) {
std::cout << "E6" << std::endl;
    return std::numeric_limits<int>::min();
  }
  if (*str == '\0') return sign*(h*3600+m*60);
  if (*str != 'm') {
std::cout << "E7" << std::endl;
    return std::numeric_limits<int>::min();
  }
  str++;
  if (*str == '\0') return sign*(h*3600+m*60);
  const int s = ParseUInt(str);
  if (s < 0) {
std::cout << "E8" << std::endl;
    return std::numeric_limits<int>::min();
  }
  if (h==90 && s > 0) {
std::cout << "E9" << std::endl;
    return std::numeric_limits<int>::min();
  }
  if (*str == '\0') {
    return sign*(h*3600+m*60+s);
  }
  if (*str != 's') {
std::cout << "E10" << std::endl;
    return std::numeric_limits<int>::min();
  }
  str++;
  if (*str == '\0') return sign*(h*3600+m*60+s);
std::cout << "E11" << std::endl;
  return std::numeric_limits<int>::min();
}

static inline
int ScanHMS(char *str) {
  return ScanDorHMS(str,'h');
}

static inline
int ScanDMS(char *str) {
  return ScanDorHMS(str,'d');
}

int main(int argc,char *argv[]) {
  int port;
  if (argc < 3 || !(std::istringstream(argv[2])>>port) ||
      port <= 0 || port >= 0x10000) {
    std::cerr << "Usage:   " << argv[0] << " <host> <port> [<ra> <dec> [<coordinate_system>]]\n"
                 "Examples: " << argv[0] << " localhost 10000 22h42m54s 56d40m06s\n"
                 "          " << argv[0] << " localhost 10000 0 0 land"
              << std::endl;
    return 1;
  }

  boost::asio::io_context io_context;

  boost::asio::ip::tcp::resolver resolver(io_context);
  boost::asio::ip::tcp::resolver::results_type endpoints =
    resolver.resolve(argv[1],argv[2]);

  boost::asio::ip::tcp::socket socket(io_context);
  boost::asio::connect(socket,endpoints);
  client
    = TcpIpConnection::Create(std::move(socket),
                              TcpIpConnection::RecvHandler(
                                &HandleStellariumTelescopeProtocolServerMsg),
                              TcpIpConnection::ErrorHandler());
  if (argc >= 5) {
    int ra = ScanHMS(argv[3]);
    int dec = ScanDMS(argv[4]);
    int coordinate_system = 0;
    if (argc >= 6) {
      coordinate_system = ParseCoordinateSystem(argv[5]);
      if (coordinate_system < 0) {
        std::cerr << "bad coordinate system \"" << argv[5] << '"'
                  << std::endl;
        return 1;
      }
    }
    if (ra != std::numeric_limits<int>::min() &&
        dec != std::numeric_limits<int>::min()) {
      const unsigned int ra_int = (((long long int)ra)<<32) / (24*3600);
      const unsigned int dec_int = (((long long int)dec)<<32) / (360*3600);
      SendMsgGoto(ra_int,dec_int,
                  static_cast<Telescope::CoordinateSystem>(coordinate_system));
    } else {
      std::cerr << "<ra> <dec> must be given like 13h2m12s -90d00m00s"
                << std::endl;
      return 1;
    }
  }

//  try {
    if (argc >= 5) {
      io_context.run_for(std::chrono::seconds(10));
    } else {
      io_context.run();
    }


//  } catch (std::exception& e) {
//    std::cerr << "Exception: " << e.what() << "\n";
//  }

  client.reset();
  return 0;
}
