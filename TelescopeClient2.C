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

#include <boost/asio.hpp>

#include <iostream>
#include <string>

static TcpIpConnection::Ptr client;

static void HandleMsgCurrentPosition(unsigned int ra_int,int dec_int) {
  std::cout << "HandleMsgCurrentPosition(" << ra_int << ',' << dec_int << ')'
            << std::endl;
}

static void SendMspGoto(unsigned int ra_int,int dec_int) {
  unsigned char msg[24];
  *reinterpret_cast<unsigned short*>(&msg) = 24;
  *reinterpret_cast<unsigned short*>(&msg[2]) = 0;
  *reinterpret_cast<unsigned long long*>(&msg[4]) = 0;
  *reinterpret_cast<unsigned int*>(&msg[12]) = ra_int;
  *reinterpret_cast<int*>(&msg[16]) = dec_int;
  *reinterpret_cast<int*>(&msg[20]) = 0;
  client->write(msg,24);
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

int main(int argc,char *argv[]) {
  int port;
  if (argc < 3 || !(std::istringstream(argv[2])>>port) ||
      port <= 0 || port >= 0x10000) {
    std::cerr << "Usage: " << argv[0] << " <host> <port> [<ra> <dec>]"
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
    double ra,dec;
    if ((std::istringstream(argv[3])>>ra) &&
        (std::istringstream(argv[4])>>dec)) {
      const unsigned int ra_int = ra * ( 2147483648.0 / 12);
      const int dec_int = dec * ( 2147483648.0 / 180);
      SendMspGoto(ra_int,dec_int);
    }
  }

//  try {
    io_context.run();
//  } catch (std::exception& e) {
//    std::cerr << "Exception: " << e.what() << "\n";
//  }

  client.reset();
  return 0;
}
