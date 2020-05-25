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

#ifndef ACCEPT_PORT_H_
#define ACCEPT_PORT_H_

#include "Connection.H"

#include <boost/asio.hpp>
#include <boost/intrusive_ptr.hpp>

#include <set>

class AcceptPort {
public:
  AcceptPort(boost::asio::io_context &io_context,unsigned short int port,
             TcpIpConnection::RecvHandler &&recv_handler)
    : acceptor(io_context),
      setup_deadline(io_context),
      port(port),
      recv_handler(std::move(recv_handler)),
      accept_connections(false) {
  }
  ~AcceptPort(void) {
  }
  void broadcast(const void *data,int size) {
    for (auto &c : connections) c->write(data,size);
  }
  void start(void) {
    if (accept_connections) return;
    accept_connections = true;
    doStart();
  }
  void stop(void) {
    if (!accept_connections) return;
    accept_connections = false;
    doStop();
  }
private:
  void doStart(void) {
    if (setupAcceptor()) {
      setupAccepting();
      return;
    }
    setup_deadline.expires_from_now(boost::posix_time::seconds(1));
    setup_deadline.async_wait(
      [this](const boost::system::error_code &e) {
        if (e) return; // timer was cancelled
        doStart();
      });
  }
  void doStop(void) {
    setup_deadline.cancel();
    boost::system::error_code ec;
    if (acceptor.is_open()) {
      acceptor.close(ec);
      if (ec) {
        std::cout << "AcceptPort(" << port << ")::stop: "
                     "close failed: "
                  << ec.message() << std::endl;
      }
    }
    for (auto &c : connections) c->hangup();
    connections.clear();
  }
  bool setupAcceptor(void) {
    if (acceptor.is_open()) return true;
    boost::system::error_code ec;
    const boost::asio::ip::tcp::endpoint ep(boost::asio::ip::tcp::v4(),port);
    acceptor.open(ep.protocol(),ec);
    if (ec) {
      std::cout << "AcceptPort(" << port << ")::setupAcceptor: "
                   "open() failed: "
                << ec.message() << std::endl;
    } else {
      acceptor.set_option(boost::asio::ip::tcp::acceptor::reuse_address(true),ec);
      if (ec) {
        std::cout << "AcceptPort(" << port << ")::setupAcceptor: "
                     "set_option(reuse_address) failed: "
                  << ec.message() << std::endl;
        acceptor.close(ec);
        return false;
      } else {
        acceptor.bind(ep,ec);
        if (ec) {
          std::cout << "AcceptPort(" << port << ")::setupAcceptor: "
                       "bind failed: "
                    << ec.message() << std::endl;
          return false;
        } else {
          acceptor.listen(3,ec);
          if (ec) {
            std::cout << "AcceptPort::stop: listen(3) failed: "
                      << ec.message() << std::endl;
          } else {
            return true;
          }
        }
      }
      acceptor.close(ec);
      if (ec) {
        std::cout << "AcceptPort(" << port << ")::setupAcceptor: "
                     "close failed: "
                  << ec.message() << std::endl;
      }
    }
    return false;
  }
  
  void setupAccepting(void) {
    acceptor.async_accept(
      [this](const boost::system::error_code& error,
             boost::asio::ip::tcp::socket &&peer) {
        if (error) {
          if (accept_connections) {
            std::cout << "AcceptPort::setupAccepting::l: accept failed: "
                      << error.message() << std:: endl;
            doStop();
            doStart();
          }
          return;
        }
        peer.set_option(boost::asio::ip::tcp::no_delay(true));
        TcpIpConnection::Ptr c(
          TcpIpConnection::Create(
            std::move(peer),
            TcpIpConnection::RecvHandler(recv_handler),
            [this](const boost::system::error_code &error,
                   TcpIpConnection &c) {
              connections.erase(TcpIpConnection::Ptr(&c));
            }));
        connections.insert(c);
        setupAccepting();
      });
  }
  int handleMsg(const unsigned char *data,int size,TcpIpConnection &c) {
      // just echo
    if (!c.write(data,size)) return -1;
    return size;
  }
private:
  boost::asio::ip::tcp::acceptor acceptor;
  boost::asio::deadline_timer setup_deadline;
  const unsigned short int port;
  std::set<TcpIpConnection::Ptr> connections;
  const TcpIpConnection::RecvHandler recv_handler;
  bool accept_connections;
};



#endif
