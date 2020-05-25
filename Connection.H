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

#ifndef CONNECTION_H_
#define CONNECTION_H_

#include "IntrusiveBase.H"

#include <boost/asio.hpp>
#include <boost/intrusive_ptr.hpp>

#include <functional>
#include <queue>
#include <iostream>

template<class T>
class Connection : public IntrusiveBase {
public:
  typedef boost::intrusive_ptr<Connection> Ptr;
    // return < 0: Error
    // return ==0: (data,size) contains no complete message,
    //             call again later
    // return > 0: a message with this amount of data has been parsed
  typedef std::function<int(const void *data,int size,
                            Connection &c)> RecvHandler;
    // called when the connection wants to die:
  typedef std::function<void(const boost::system::error_code &error,
                             Connection &c) > ErrorHandler;
  static Ptr Create(T &&peer,
                    RecvHandler &&recv_handler,
                    ErrorHandler &&error_handler) {
    return new Connection(std::move(peer),
                          std::move(recv_handler),
                          std::move(error_handler));
  }
  bool write(const void *data,int size) {
    if (size > (int)sizeof(send_buf) - send_buf_end) return false;
    memcpy(send_buf+send_buf_end,data,size);
    if (send_buf_end == 0) {
      send_buf_written = 0;
      send_buf_end = size;
      callAsyncWrite();
    } else {
      send_buf_end += size;
    }
    return true;
  }
//  bool write(const std::string &data) {
//    return write(data.data(),data.size());
//  }
  void hangup(void) {
    boost::system::error_code ec;
    peer.close(ec);
    if (ec) {
      std::cout << "Connection::hangup: close failed: " << ec.message()
                << std::endl;
    }
  }
protected:
  Connection(T &&peer,
             RecvHandler &&recv_handler,
             ErrorHandler &&error_handler)
    : peer(std::move(peer)),
      recv_handler(std::move(recv_handler)),
      error_handler(std::move(error_handler)) {
//    std::cout << "Connection<" << typeid(T).name() << ">::Connection" << std::endl;
    read();
  }
  ~Connection(void) {
//    std::cout << "Connection<" << typeid(T).name() << ">::~Connection" << std::endl;
  }
  void resetBuffers(void) {
    recv_used = 0;
    send_buf_end = 0;
    send_buf_written = 0;
  }
  void callAsyncWrite(void) {
    peer.async_write_some(
      boost::asio::buffer(send_buf+send_buf_written,
                          send_buf_end-send_buf_written),
      [self=Ptr(this)](const boost::system::error_code &error,
                       std::size_t size) {
        if (error) {
          std::cout << "Connection::write::l: write failed: "
                    << error.message() << std:: endl;
          self->hangup(); // close socket, cancel io
          self->error_handler(error,*self);
            // do not write again, the connection shall die.
        } else {
          self->send_buf_written += size;
          if (self->send_buf_written >= self->send_buf_end) {
            self->send_buf_end = 0;
          } else {
            if (self->send_buf_written > (int)sizeof(self->send_buf) / 2) {
                // rearrange: make space
              self->send_buf_end -= self->send_buf_written;
              memmove(self->send_buf,
                      self->send_buf+self->send_buf_written,
                      self->send_buf_end);
              self->send_buf_written = 0;
            }
            self->callAsyncWrite();
          }
        }
      });
  }
  void read(void) {
    peer.async_read_some(
      boost::asio::buffer(recv_buf+recv_used,sizeof(recv_buf)-recv_used),
      [self=Ptr(this)](const boost::system::error_code &error,
                       std::size_t size) {
        if (error) {
          std::cout << "Connection::read::l: read failed: "
                    << error.message() << std:: endl;
          self->hangup(); // close socket, cancel io
          self->error_handler(error,*self);
            // do not read again, the connection shall die.
        } else {
          self->recv_used += size;
          int start = 0;
          for (;;) {
            const int rc = self->recv_handler(self->recv_buf+start,
                                              self->recv_used-start,
                                              *self);
            if (rc == 0) {
              self->recv_used -= start;
              memmove(self->recv_buf,
                      self->recv_buf+start,
                      self->recv_used);
              break;
            } else if (0 < rc) {
              start += rc;
              if (start >= self->recv_used) {
                self->recv_used = 0;
                break;
              }
            } else {
              std::cout << "Connection::read::l: handleMsg failed: "
                        << rc << std::endl;
                // do not read again
              self->hangup(); // close socket, cancel io
              self->error_handler(error,*self);
            }
          }
          self->read();
        }
      });
  }
private:
  T peer;
  unsigned char recv_buf[1024];
  unsigned char send_buf[1024];
  int recv_used = 0;
  int send_buf_end = 0;
  int send_buf_written = 0;
  const RecvHandler recv_handler;
  const ErrorHandler error_handler;
};


/*
RquRspConnection is actually not good enough for serial communication,
I do not need it.

template<class T>
class RquRspConnection : public Connection<T> {
public:
  typedef boost::intrusive_ptr<RquRspConnection> Ptr;
  static Ptr Create(T &&peer) {
    return new RquRspConnection(std::move(peer));
  }
  typedef std::function<int(const void *data,int size)> ResponseFunc;
  bool request(const void *data,int size,ResponseFunc &&f) {
    if (!write(data,size)) return false;
    response_queue.emplace(std::move(f));
    return true;
  }
//  bool request(const std::string &data,ResponseFunc &&f) {
//    return request(data.data(),data.size(),std::move(f));
//  }
  using Connection<T>::write;
private:
  RquRspConnection(T &&peer)
    : Connection<T>(
        std::move(peer),
        [this](const void *data,int size,
               Connection<T> &c)->int {
          const int rc = response_queue.front()(data,size);
          if (rc > 0) {
            response_queue.pop();
            return rc;
          }
          if (rc < 0) {
            do {
              response_queue.front()(0,0);
              response_queue.pop();
            } while (!response_queue.empty());
            Connection<T>::resetBuffers();
          }
          return rc;
        },
        [this](const boost::system::error_code &error,
               Connection<T> &c) {
          std::cout << "RquRspConnection: error "
                    << error.message() << std::endl;
          throw nullptr;
        }) {}
  std::queue<ResponseFunc> response_queue;
};
*/

typedef Connection<boost::asio::ip::tcp::socket> TcpIpConnection;
//typedef Connection<boost::asio::serial_port> SerialConnection;
//typedef RquRspConnection<boost::asio::ip::tcp::socket> TcpIpRquRspConnection;
//typedef RquRspConnection<boost::asio::serial_port> SerialRquRspConnection;

#endif
