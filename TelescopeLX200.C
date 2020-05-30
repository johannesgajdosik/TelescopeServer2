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
#include "Connection.H"
#include "PrintRaDec.H"
#include "Time.H"

#include "Vector.H"

#include <boost/asio.hpp>

#include <iostream>
#include <iomanip>


#include <time.h>
#ifdef _WIN32
  #define TZSET _tzset
  #define TZNAME _tzname
  #define TIMEZONE _timezone
  #define DAYLIGHT _daylight
#else
  #define TZSET tzset
  #define TZNAME tzname
  #define TIMEZONE timezone
  #define DAYLIGHT daylight
#endif




struct LX200RaDec {
  LX200RaDec(void) {}
  LX200RaDec(unsigned int ra,bool dec_sign,unsigned int dec) // seconds
    : ra(ra),dec_sign(dec_sign),dec(dec) {
    if (dec > 90*3600) abort();
  }
  unsigned int getRa(void) const {return ra;}
  bool getDecSign(void) const {return dec_sign;}
  unsigned int getAbsDec(void) const {return dec;}
  int getDec(void) const {return (dec_sign ? -dec : dec);}
  LX200RaDec(unsigned int ra_int,int dec_int) {
    set(ra_int,dec_int);
  }
  void set(unsigned int ra_int,int dec_int) {
    ra  = ((unsigned int)((                                                      ra_int  * (24*3600ULL) + 0x80000000ULL) >> 32));
    dec_sign = (dec_int < 0);
    dec = ((unsigned int)((  (dec_sign ? -(unsigned int)dec_int : (unsigned int)dec_int) * (360*3600ULL) + 0x80000000ULL) >> 32));
  }
  unsigned int getRaInt(void) const {
    return ((((unsigned long long int)ra)<<32)+(12*3600ULL))/(24*3600ULL);
  }
  double getRaRad(void) const {return ra*(M_PI/(12*3600.0));}
  int getDecInt(void) const {
    int dec_int = ((((unsigned long long int)dec)<<32)+(180*3600ULL))/(360*3600ULL);
    if (dec_sign) dec_int = -dec_int;
    return dec_int;
  }
  double getDecRad(void) const {
    const double dec_rad = dec*(M_PI/(180*3600.0));
    return (dec_sign ? -dec_rad : dec_rad);
  }
  unsigned int ra;
  bool dec_sign;
  unsigned int dec;
};

static
std::ostream &operator<<(std::ostream &o,const LX200RaDec &x) {
  o << "(Ra:" << PrintRaSeconds(x.getRa())
    << ",Dec:" << PrintDecSeconds(x.getDecSign(),x.getAbsDec()) << ')';
  return o;
}









class TelescopeLX200 : public Telescope {
public:
  static Ptr Create(const std::string &args,
                    OpenedClosedFunction &&opened_closed,
                    PositionFunction &&announce_position,
                    boost::asio::io_context &io_context);
  class Command {
  public:
    virtual ~Command(void) {}
    virtual void execAsync(void) = 0;
    virtual unsigned int getTimeoutMicros(void) const {return 500000;}
    virtual void print(std::ostream &o) const = 0;
  protected:
    Command(TelescopeLX200 &telescope) : telescope(telescope) {}
    TelescopeLX200 &telescope;
  };
  TelescopeLX200(const TelescopeLX200&) = delete;
  TelescopeLX200 &operator=(const TelescopeLX200&) = delete;
private:
  TelescopeLX200(const std::string &args,
                   OpenedClosedFunction &&opened_closed,
                   PositionFunction &&announce_position,
                   boost::asio::io_context &io_context);
  ~TelescopeLX200(void);
  bool initializationOk(void) const override {return true;}

    // RS232 send/receive
  void recvRsp(std::function<int(void)> &&rsp_data_received);
  void sendRqu(const char *data,int size,
               std::function<int(void)>  &&rsp_data_received);
  void sendMsg(const char *data,int size,
               std::function<void(void)> &&finished);

    // initialization
  void initialize(void);
  void initPrecessionMatrix(void);
  class CommandInit;

    // geographic location
  void locationReceived(unsigned int longitude_int,int latitude_int);
  class CommandGetLoc;
  void getLoc(void);

    // telescope pointing position in the sky
  void positionReceived(const LX200RaDec &ra_dec);
  class CommandGetPos;
  void getPos(void);

    // Goto
  class CommandGoto;
  void gotoPosition(const LX200RaDec &ra_dec);
  void gotoPosition(const unsigned int ra_int_j2000,const int dec_int_j2000) override;

    // move
  class CommandMove;
  void move(short int horz,short int vert,
            unsigned int validity_micros) override;
  
    // guide
  class CommandGuide;
  void guide(int d_ra_micros,int d_dec_micros) override;

    // command handling
  void init(unsigned int drain_micros = 2000000);
  void commandFinished(void);
  void doSomething(void);
private:
  const std::string args;
  boost::asio::serial_port serial;
  boost::asio::deadline_timer serial_deadline;
  boost::asio::deadline_timer command_deadline;
  boost::asio::deadline_timer get_pos_deadline;
  boost::asio::deadline_timer move_deadline;
  std::unique_ptr<Command> curr_command;
  std::unique_ptr<CommandInit> next_command_init;
  std::unique_ptr<CommandGetLoc> next_command_get_loc;
  std::unique_ptr<CommandGetPos> next_command_get_pos;
  std::unique_ptr<CommandGoto> next_command_goto;
  std::unique_ptr<CommandMove> next_command_move;
  std::unique_ptr<CommandGuide> next_command_guide;

  Matrix<double,3,3> precession_matrix;
  Matrix<double,3,3> geographic_pos_orientation;
  unsigned int longitude_int = 0;
  int latitude_int = 0x7FFFFFFF;
  long long int end_of_last_goto = (-0x7FFFFFFFFFFFFFFFLL-1LL);

    // max buffer contents:
    // 0         1         2         3         4         5         6         
    // 0123456789012345678901234567890123456789012345678901234567890123456789
    // 1Updating        planetary data. #                                #
  char recv_buf[128];
  int recv_used;
  short int curr_horz = 0;
  short int curr_vert = 0;
  short int curr_speed = 0;

  bool telescope_online = false;
  bool language_v3 = false;
};

static inline
std::ostream &operator<<(std::ostream &o,const TelescopeLX200::Command &c) {
  c.print(o);
  return o;
}

static const bool telescope_LX200_registered
  = Telescope::RegisterCreationFunction("LX200",TelescopeLX200::Create);




//------------------------------------------------------------------------------





//-----construction/destruction---------------------------------------------------


Telescope::Ptr
TelescopeLX200::Create(const std::string &args,
                         OpenedClosedFunction &&opened_closed,
                         PositionFunction &&announce_position,
                         boost::asio::io_context &io_context) {
  return new TelescopeLX200(args,
                              std::move(opened_closed),
                              std::move(announce_position),
                              io_context);
}

TelescopeLX200::TelescopeLX200(const std::string &args,
                                   OpenedClosedFunction &&opened_closed,
                                   PositionFunction &&announce_position,
                                   boost::asio::io_context &io_context)
                 :Telescope(std::move(opened_closed),std::move(announce_position)),
                  args(args),
                  serial(io_context),
                  serial_deadline(io_context),
                  command_deadline(io_context),
                  get_pos_deadline(io_context),
                  move_deadline(io_context),
                  language_v3(false) {
  io_context.dispatch(std::bind(&TelescopeLX200::initialize,this));
  std::cout << PrintTime() << " "
               "TelescopeLX200::TelescopeLX200(" << args << ')' << std::endl;
}

TelescopeLX200::~TelescopeLX200(void) {
  if (telescope_online) {
    opened_closed(false);
    telescope_online = false;
  }
  std::cout << PrintTime() << " "
               "TelescopeLX200::~TelescopeLX200" << std::endl;
}

//------------------------------------------------------------------------------




//-------------------RS232 communication------------------------------------


void TelescopeLX200::recvRsp(std::function<int(void)> &&rsp_data_received) {
  const int to_read = (int)(sizeof(recv_buf)-1) - recv_used;
  if (to_read <= 0) {
    std::cout << PrintTime() << " "
                 "TelescopeLX200::recvRsp: recv_buf full, discarding response"
              << std:: endl;
    init();
    return;
  }
  serial.async_read_some(
    boost::asio::buffer(recv_buf+recv_used,to_read),
    [this,f = std::move(rsp_data_received)](
        const boost::system::error_code &error,std::size_t size) mutable {
      if (error) {
        if (error == boost::system::errc::io_error) {
          std::cout << PrintTime() << " "
                       "TelescopeLX200::recvRsp::l: "
                       "IO ERROR" << std::endl;
          initialize();
          return;
        }
        if (error == boost::asio::error::eof) {
          std::cout << PrintTime() << " "
                       "TelescopeLX200::recvRsp::l: "
                       "EOF" << std::endl;
          initialize();
          return;
        }
        std::cout << PrintTime() << " "
                     "TelescopeLX200::recvRsp::l: read failed: "
                  << error.message() << std:: endl;
        init();
        return;
      }
      recv_used += size;
      recv_buf[recv_used] = '\0';
      const int rc = f();
      if (rc == 0) {
        std::function<int(void)> f2;
        f2.swap(f);
        recvRsp(std::move(f2));
        return;
      }
        // f() has alread finished the command
      if (rc < 0) {
        init();
      }
    });
}

void TelescopeLX200::sendRqu(const char *data,int size,
                             std::function<int(void)>  &&rsp_data_received) {
  sendMsg(data,size,
          [this,f = std::move(rsp_data_received)]() mutable {
            recv_used = 0;
            std::function<int(void)> f2;
            f2.swap(f);
            recvRsp(std::move(f2));
          });
}

void TelescopeLX200::sendMsg(const char *data,int size,
                             std::function<void(void)> &&finished) {
  serial.async_write_some(
    boost::asio::buffer(data,size),
    [this,data,size,f = std::move(finished)] (
        const boost::system::error_code &error,std::size_t written) mutable {
      if (error) {
        if (error == boost::system::errc::io_error) {
          std::cout << PrintTime() << " "
                       "TelescopeLX200::sendMsg::l: "
                       "IO ERROR" << std::endl;
          initialize();
          return;
        }
        std::cout << PrintTime() << " "
                     "TelescopeLX200::sendMsg::l: write failed: "
                  << error.message() << std:: endl;
        init();
        return;
      }
      data += written;
      size -= written;
      if (size > 0) {
        sendMsg(data,size,std::move(f));
        return;
      }
      f();
    });
}

//------------------------------------------------------------------------------





//------------------Initialization--------------------------------------

static
void Lea406aPrecMat(const double t,double prec[9]) {
  const double prec_dz =
    (((((-0.0002*t-0.0285)*t-0.0583)*t+18.0183)*t +30.2226)*t+23060.9097)*t
      *(M_PI/(180*3600));
  const double cos_dz = cos(prec_dz);
  const double sin_dz = sin(prec_dz);
  const double prec_z =
    (((((-0.0001*t-0.0301)*t-0.2821)*t+18.2667)*t+109.5270)*t+23060.9097)*t
      *(M_PI/(180*3600));
  const double cos_z = cos(prec_z);
  const double sin_z = sin(prec_z);
  const double prec_t =
    ((((( 0.0004*t-0.0127)*t-0.0731)*t-41.8238)*t -42.6566)*t+20042.0207)*t
      *(M_PI/(180*3600));
  const double cos_t = cos(prec_t);
  const double sin_t = sin(prec_t);
  prec[0] =  -sin_dz * sin_z + cos_dz * cos_z * cos_t;
  prec[1] =   sin_dz * cos_z + cos_dz * sin_z * cos_t;
  prec[2] =                    cos_dz         * sin_t;
  prec[3] =  -cos_dz * sin_z - sin_dz * cos_z * cos_t;
  prec[4] =   cos_dz * cos_z - sin_dz * sin_z * cos_t;
  prec[5] =                   -sin_dz         * sin_t;
  prec[6] =                            -cos_z * sin_t;
  prec[7] =                            -sin_z * sin_t;
  prec[8] =                                     cos_t;
}

void TelescopeLX200::initPrecessionMatrix(void) {
  const double jd_minus_j2000 = TimeToJD_minus_J2000(GetNow());
  const double t = jd_minus_j2000 * (1.0 / 365250.0);
  Lea406aPrecMat(t,precession_matrix);
}

void TelescopeLX200::initialize(void) {
  if (telescope_online) {
    telescope_online = false;
    opened_closed(false);
  }
  next_command_init.reset();
  next_command_get_loc.reset();
  next_command_get_pos.reset();
  next_command_goto.reset();
  next_command_move.reset();
  next_command_guide.reset();

  boost::system::error_code ec;
  serial.close(ec);
  if (ec) {
    std::cout << PrintTime() << " "
                 "TelescopeLX200::initialize: serial port::close failed: "
              << ec.message() << std::endl;
  }
  serial_deadline.cancel();
  command_deadline.cancel();
  get_pos_deadline.cancel();
  move_deadline.cancel();
  curr_command.reset();

  serial.open(args,ec);
  if (ec) {
    std::cout << PrintTime() << " "
                 "TelescopeLX200::initialize: serial port::open(" << args << ") failed: "
              << ec.message() << std::endl;
    serial_deadline.expires_from_now(boost::posix_time::microseconds(2000000));
    serial_deadline.async_wait(
      [this](const boost::system::error_code &e) {
        if (e) return; // timer was cancelled
        initialize();
      });
    return;
  }
  serial.set_option(boost::asio::serial_port_base::baud_rate(9600));
  serial.set_option(boost::asio::serial_port_base::character_size(8));
  serial.set_option(boost::asio::serial_port_base::stop_bits(boost::asio::serial_port_base::stop_bits::one));
  serial.set_option(boost::asio::serial_port_base::parity(boost::asio::serial_port_base::parity::none));
  serial.set_option(boost::asio::serial_port_base::flow_control(boost::asio::serial_port_base::flow_control::none));
  std::cout << PrintTime() << " "
               "TelescopeLX200::initialize: serial port " << args << "opened" << std::endl;
  initPrecessionMatrix();
  init(0);
}


class TelescopeLX200::CommandInit : public TelescopeLX200::Command {
public:
  CommandInit(TelescopeLX200 &telescope,unsigned int drain_micros)
    : Command(telescope),drain_deadline(telescope.serial.get_executor()),
      drain_micros(drain_micros) {}
  void execAsync(void) override {
    sendStopRqu();
  }
  void print(std::ostream &o) const override {
    o << "CommandInit(" << drain_micros << ')';
  }
private:
  unsigned int getTimeoutMicros(void) const override {return 10*1000000 + drain_micros;}
  void sendStopRqu(void) {
    if (telescope.telescope_online) {
      telescope.telescope_online = false;
      telescope.opened_closed(false);
    }
    if (drain_micros == 0) {
      sendGetCalendarFormatRqu();
      return;
    }
    drain_deadline.expires_from_now(boost::posix_time::microseconds(drain_micros));
    drain_deadline.async_wait(
      [t = &telescope](const boost::system::error_code &error) {
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandInit::sendStopRqu::l: "
                     "draining deadline reached, " << error.message()
                  << std:: endl;
          // cancel serial read/write requests regardless of error:
        boost::system::error_code ec;
        t->serial.cancel(ec);
        if (ec) {
          std::cout << PrintTime() << " "
                       "TelescopeLX200::CommandInit::sendStopRqu::l: "
                       "serial.cancel failed: " << ec.message()
                    << std:: endl;
        }
      });
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandInit::sendStopRqu: "
                 "sending :Q#, draining for " << drain_micros << "us" << std::endl;
    telescope.sendMsg(":Q#",3,std::bind(&TelescopeLX200::CommandInit::recvStopRsp,this));
  }
  void recvStopRsp(void) {
    telescope.serial.async_read_some(
      boost::asio::buffer(telescope.recv_buf,sizeof(telescope.recv_buf)-1),
      [this](const boost::system::error_code &error,std::size_t size) {
        if (error) {
          if (error == boost::system::errc::io_error) {
            std::cout << PrintTime() << " "
                         "TelescopeLX200::CommandInit::recvStopRsp::l: "
                         "IO ERROR" << std::endl;
            telescope.initialize();
            return;
          }
          if (error == boost::asio::error::eof) {
            std::cout << PrintTime() << " "
                         "TelescopeLX200::CommandInit::recvStopRsp::l: "
                         "EOF" << std::endl;
            telescope.initialize();
            return;
          }
          if (error == boost::asio::error::operation_aborted) {
            std::cout << PrintTime() << " "
                         "TelescopeLX200::CommandInit::recvStopRsp::l: "
                         "draining finished"
                      << std:: endl;
            telescope.recv_used = 0;
            sendGetCalendarFormatRqu();
            return;
          }
          std::cout << PrintTime() << " "
                       "TelescopeLX200::CommandInit::recvStopRsp::l: "
                       "unexpected error: " << error.message() << std::endl;
          telescope.initialize();
          return;
          
        }
        telescope.recv_buf[size] = '\0';
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandInit::recvStopRsp::l: "
                     "received \"" << telescope.recv_buf << "\", draining..."
                  << std::endl;
        recvStopRsp();
      });
  }
  void sendGetCalendarFormatRqu(void) {
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandInit::sendGetCalendarFormatRqu: "
                 "sending :Gc#" << std::endl;
    telescope.sendRqu(":Gc#",4,std::bind(&TelescopeLX200::CommandInit::recvGetCalendarFormatRsp,this));
  }
  int recvGetCalendarFormatRsp(void) {
    if (telescope.recv_used < 3) return 0;
    const char *p = telescope.recv_buf;
    const bool expect_bracket = (*p == '(');
    if (expect_bracket) {
      p++;
      if (telescope.recv_used < 5) return 0;
    }
    const char *const format = p;
    p += 2;
    if (expect_bracket && *p++ != ')') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandInit::recvGetCalendarFormatRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                   "error: ')' expected" << std::endl;
      return -1;
    }
    if (*p++ != '#') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandInit::recvGetCalendarFormatRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                   "error: '#' expected" << std::endl;
      return -1;
    }
    if (0 == strncmp(format,"12",2)) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandInit::recvGetCalendarFormatRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                   "12h-Format received" << std::endl;
      sendToggleTimeFormatMsg();
      return (p-telescope.recv_buf);
    }
    if (0 == strncmp(format,"24",2)) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandInit::recvGetCalendarFormatRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                   "24h-Format received" << std::endl;
      sendSetTimeZoneRqu();
      return (p-telescope.recv_buf);
    }
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandInit::recvGetCalendarFormatRsp: \""
              << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                 "error: bad format received" << std::endl;
    return -1;
  }
  void sendToggleTimeFormatMsg(void) {
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandInit::sendToggleTimeFormatMsg: "
                 "sending :H#" << std::endl;
    telescope.sendMsg(":H#",3,std::bind(&TelescopeLX200::CommandInit::sendSetTimeZoneRqu,this));
  }
  void sendSetTimeZoneRqu(void) {
    TZSET();
    int time_zone = TIMEZONE;
    int is_dst = DAYLIGHT;
      // :SGsHH.H#
      // Set the number of hours added to local time to yield UTC
    gmt_offset = time_zone - (is_dst ? 3600 : 0);
    buf[0] = ':';
    buf[1] = 'S';
    buf[2] = 'G';
    int x = gmt_offset;
    if (x < 0) {buf[3] = '-';x = -x;}
    else {buf[3] = '+';}
    x /= 60;
    if (x > 24*60) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandInit::sendSetTimeZoneRqu: "
                   "bad gmt_offset: " << gmt_offset << ". "
                   "Assuming 0. Fix the code." << std::endl;
      x = 0;
      buf[3] = '+';
    }
    x /= 6;
    int y;
    y = x / 10; buf[7] = '0' + (char)(x-10*y);
    buf[6] = '.';
    x = y / 10; buf[5] = '0' + (char)(y-10*x);
                buf[4] = '0' + (char)(x);
    buf[8] = '#';
    buf[9] = '\0';
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandInit::sendSetTimeZoneRqu: "
                 "sending " << buf << std::endl;
    telescope.sendRqu(buf,9,std::bind(&TelescopeLX200::CommandInit::recvSetTimeZoneRsp,this));
  }
  int recvSetTimeZoneRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandInit::recvSetTimeZoneRsp \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                   "Bad Response" << std::endl;
      return -1;
    }
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandInit::recvSetTimeZoneRsp \""
              << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                 "ok" << std::endl;
    sendSetTimeRqu();
    return 1;
  }
  void sendSetTimeRqu(void) {
    unsigned long long int x = GetNow() / 1000000LL - gmt_offset;
    buf[0] = ':';
    buf[1] = 'S';
    buf[2] = 'L';
    unsigned long long int y;
    y = x / 10; buf[10] = '0' + (char)(x-10ULL*y);
    x = y /  6; buf[ 9] = '0' + (char)(y- 6ULL*x);
    buf[8] = ':';
    y = x / 10; buf[ 7] = '0' + (char)(x-10ULL*y);
    x = y /  6; buf[ 6] = '0' + (char)(y- 6ULL*x);
    buf[5] = ':';
    days = x / 24; x -= 24*days;
    y = x / 10; buf[ 4] = '0' + (char)(x-10ULL*y);
                buf[ 3] = '0' + (char)(        y);
    buf[11] = '#';
    buf[12] = '\0';
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandInit::sendSetTimeRqu: "
                 "sending " << buf << std::endl;
    telescope.sendRqu(buf,12,std::bind(&TelescopeLX200::CommandInit::recvSetTimeRsp,this));
  }
  int recvSetTimeRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandInit::recvSetTimeRsp \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                   "Bad Response" << std::endl;
      return -1;
    }
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandInit::recvSetTimeRsp \""
              << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                 "ok" << std::endl;
    sendSetDateRqu();
    return 1;
  }
  void sendSetDateRqu(void) {
    int y,m,d;
    DecomposeDay(days,y,m,d);
    m++;d++;
    int x;
    buf[0] = ':';
    buf[1] = 'S';
    buf[2] = 'C';
    x = m / 10; buf[ 4] = '0' + (char)(m-10*x);
                buf[ 3] = '0' + (char)(     x);
    buf[5] = '/';
    x = d / 10; buf[ 7] = '0' + (char)(d-10*x);
                buf[ 6] = '0' + (char)(     x);
    buf[8] = '/';
    x = y / 10; buf[10] = '0' + (char)(y-10*x);
    y = x / 10; buf[ 9] = '0' + (char)(x-10*y);
    buf[11] = '#';
    buf[12] = '\0';
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandInit::sendSetDateRqu: "
                 "sending " << buf << std::endl;
    telescope.sendRqu(buf,12,std::bind(&TelescopeLX200::CommandInit::recvSetDateRsp,this));
  }
  int recvSetDateRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandInit::recvSetDateRsp \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                   "Bad Response" << std::endl;
      return -1;
    }
    int count = 0;
    for (int i=1;i<telescope.recv_used;i++) {
      if (telescope.recv_buf[i] == '#') {
        count++;
        if (count == 2) {
          std::cout << PrintTime() << " "
                       "TelescopeLX200::CommandInit::recvSetDateRsp ok: \""
                    << std::string(telescope.recv_buf,telescope.recv_used) << '"' << std::endl;
          telescope.getLoc();
          if (!telescope.telescope_online) {
            telescope.telescope_online = true;
            telescope.opened_closed(true);
          }
          telescope.commandFinished();
          return (i+1);
        }
      }
    }
    return 0;
  }
private:
  boost::asio::deadline_timer drain_deadline;
  const unsigned int drain_micros;
  int gmt_offset,days;
  char buf[13];
};


//------------------------------------------------------------------------




static void FillDecimal(char *data,unsigned int size,unsigned int x) {
  data += size;
  while (size > 0) {
    const unsigned int h = x/10;
    *--data = '0' + (x-10*h);
    x = h;
    size--;
  }
}

static bool ParseDecimal(const char *data,unsigned int size,unsigned int &x) {
  x = 0;
  while (size > 0) {
    const char c = *data++;
    if (c < '0' || '9' < c) return false;
    x = 10u*x + (unsigned int)(c-'0');
    size--;
  }
  return true;
}



//----------Geographic Location----------------------------------------

static
void SetLongLat(double longitude,double latitude,
                Matrix<double,3,3> &orientation) {
  const double cl = cos(longitude);
  const double sl = sin(longitude);
  const double cp = cos(latitude);
  const double sp = sin(latitude);
  orientation =
    Matrix<double,3,3>(
      sp, 0,-cp,
       0, 1,  0,
      cp, 0, sp)
  * Matrix<double,3,3>(
       sl, cl, 0,
      -cl, sl, 0,
        0,  0, 1);
}

void TelescopeLX200::locationReceived(unsigned int new_longitude_int,int new_latitude_int) {
  if (longitude_int != new_longitude_int ||
      latitude_int != new_latitude_int) {
    longitude_int = new_longitude_int;
    latitude_int = new_latitude_int;
    std::cout << PrintTime() << " "
                 "TelescopeLX200::locationReceived: new location ("
              << PrintRaInt(longitude_int) << ','
              << PrintDecInt(latitude_int) << ')' << std::endl;
    SetLongLat(longitude_int*(M_PI/2147483648.0),
               latitude_int*(M_PI/2147483648.0),
               geographic_pos_orientation);
  }
}

class TelescopeLX200::CommandGetLoc : public TelescopeLX200::Command {
public:
  CommandGetLoc(TelescopeLX200 &telescope) : Command(telescope) {}
  void execAsync(void) override {
    sendGetLatitudeRqu();
  }
  void print(std::ostream &o) const override {o << "CommandGetLoc";}
private:
  void sendGetLatitudeRqu(void) {
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandInit::sendGetLatitudeRqu: "
                 "sending :Gt#" << std::endl;
    telescope.sendRqu(":Gt#",4,std::bind(&TelescopeLX200::CommandGetLoc::recvGetLatitudeRsp,this));
  }
  int recvGetLatitudeRsp(void) {
    if (telescope.recv_used < 7) return 0;
    if (telescope.recv_buf[6] != '#') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLatitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": no '#' at end, ignoring response" << std::endl;
      return -1;
    }
    const bool latitude_sign = (telescope.recv_buf[0] == '-');
    if (!latitude_sign && (telescope.recv_buf[0] != '+')) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLatitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": latitude sign must be '+' or '-', ignoring response" << std::endl;
      return -1;
    }
    unsigned int degrees;
    if (!ParseDecimal(telescope.recv_buf+1,2,degrees)) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLatitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad degrees: no number, ignoring response" << std::endl;
      return -1;
    }
      // spec: '*', scope: 223
    if (telescope.recv_buf[3] != (char)223 && telescope.recv_buf[3] != '*') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLatitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": no degree symbol between degrees and minutes, ignoring response" << std::endl;
      return -1;
    }
    unsigned int minutes;
    if (!ParseDecimal(telescope.recv_buf+4,2,minutes)) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLatitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad minutes: no number, ignoring response" << std::endl;
      return -1;
    }
    if (minutes >= 60) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLatitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad minutes: too big, ignoring response" << std::endl;
      return -1;
    }
    minutes += 60*degrees;
    if (minutes > 60*90) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLatitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad degrees: too big, ignoring response" << std::endl;
      return -1;
    }
    latitude_int = ((((unsigned long long int)minutes)<<32)+180*60ULL) / (360*60ULL);
    if (latitude_sign) latitude_int = -latitude_int;
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandGetLoc::recvGetLatitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": ok: " << PrintDecInt(latitude_int) << std::endl;
    sendGetLongitudeRqu();
    return 7;
  }
  void sendGetLongitudeRqu(void) {
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandInit::sendGetLongitudeRqu: "
                 "sending :Gg#" << std::endl;
    telescope.sendRqu(":Gg#",4,std::bind(&TelescopeLX200::CommandGetLoc::recvGetLongitudeRsp,this));
  }
  int recvGetLongitudeRsp(void) {
    if (telescope.recv_used < 7) return 0;
      // the new spec says there is a sign, but my scope sends no sign.
    const bool longitude_sign = (telescope.recv_buf[0] == '-');
    const char *p = telescope.recv_buf;
    if (longitude_sign || telescope.recv_buf[0] == '+') {
      if (telescope.recv_used < 8) return 0;
      p++;
    }
    if (p[6] != '#') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLongitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": no '#' at end, ignoring response" << std::endl;
      return -1;
    }
    unsigned int degrees;
    if (!ParseDecimal(p,3,degrees)) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLongitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad degrees: no number, ignoring response" << std::endl;
      return -1;
    }
    if (degrees >= 360) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLongitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad degrees: too big, ignoring response" << std::endl;
      return -1;
    }
      // spec: '*', scope: 223
    if (p[3] != (char)223 && p[3] != '*') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLongitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": no degree symbol between degrees and minutes, ignoring response" << std::endl;
      return -1;
    }
    unsigned int minutes;
    if (!ParseDecimal(p+4,2,minutes)) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLongitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad minutes: no number, ignoring response" << std::endl;
      return -1;
    }
    if (minutes >= 60) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetLoc::recvGetLongitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad minutes: too big, ignoring response" << std::endl;
      return -1;
    }
    unsigned int longitude_int = ((((unsigned long long int)(degrees*60+minutes))<<32)+180*60ULL) / (360*60ULL);
    if (longitude_sign) longitude_int = -longitude_int;
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandGetLoc::recvGetLongitudeRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": ok: " << PrintRaInt(longitude_int) << std::endl;
    telescope.locationReceived(longitude_int,latitude_int);
    telescope.getPos();
    telescope.commandFinished();
    return (p-telescope.recv_buf);
  }
  int latitude_int;
};

void TelescopeLX200::getLoc(void) {
  if (!next_command_get_loc) {
    next_command_get_loc = std::make_unique<CommandGetLoc>(*this);
  }
  doSomething();
}


//------------------------------------------------------------------------



//---------------Pointing Position----------------------------------------

void TelescopeLX200::positionReceived(const LX200RaDec &ra_dec) {
  const double ra0 = ra_dec.getRaRad();
  const Vector<double,3> v0 = PolarToRect(ra0,ra_dec.getDecRad());
  const Vector<double,3> v = precession_matrix*v0;
  double ra_j2000,dec_j2000;
  RectToPolar(v,ra_j2000,dec_j2000);
  const unsigned int  ra_int_j2000 = (unsigned int)floor(0.5+ ra_j2000*(2147483648.0/M_PI));
  const          int dec_int_j2000 =          (int)floor(0.5+dec_j2000*(2147483648.0/M_PI));
  announce_position(ra_int_j2000,dec_int_j2000);

    // from here on: check if telescope might bumps into gate
  const long long int now = GetNow();
  if (now < end_of_last_goto) return;
  const double jd_minus_j2000 = TimeToJD_minus_J2000(now);
  constexpr double W0 = (190.16+90)/360.0;
  constexpr double dW = 360.9856235/360.0;
  double W = W0 + dW * jd_minus_j2000;
  W -= floor(W);
  const double W_rad = W* (2.0*M_PI);
  const double cW = cos(W_rad);
  const double sW = sin(W_rad);
  Matrix<double,3,3> earth_orientation
    = Matrix<double,3,3>( sW,-cW,0.0,
                          cW, sW,0.0,
                         0.0,0.0,1.0);
  if (latitude_int > 0x40000000) {
    std::cout << PrintTime() << " "
                 "TelescopeLX200::positionReceived" << ra_dec << ": "
                 "J2000(Ra:" << PrintRaInt(ra_int_j2000)
              << ",Dec:" << PrintDecInt(dec_int_j2000) << "): "
                 "geographic location unknown, cannot check for forbidden position"
              << std::endl;
  } else {
    const unsigned int W_int = (unsigned int)((unsigned long long int)(W*4294967296.0+0.5));
    const unsigned int hour_angle_int = W_int-longitude_int-ra_dec.getRaInt();

    Matrix<double,3,3> equat_orientation = geographic_pos_orientation
                                         * earth_orientation;
    double az,alt;
    RectToPolar(equat_orientation*v0,az,alt);
#define PRINT_PERIODIC_POSITION

#ifdef PRINT_PERIODIC_POSITION
      std::cout << PrintTime() << " "
                   "TelescopeLX200::positionReceived" << ra_dec << ": "
                   "J2000(Ra:" << PrintRaInt(ra_int_j2000)
                << ",Dec:" << PrintDecInt(dec_int_j2000) << "), "
                   "H:" << PrintRaInt(hour_angle_int)
                << ", Az:" << PrintRaRad(az) << ", Alt:" << PrintDecRad(alt);
#endif

    if (alt < -2.0*(M_PI/180)) {
      alt = 60*(M_PI/180.0);
      az = (360-90)*(M_PI/180.0);
      double ra,dec;
      RectToPolar(PolarToRect(az,alt)*equat_orientation,ra,dec);
      const unsigned int goto_ra_int = (unsigned int)floor(0.5+ra*(2147483648.0/M_PI));
      const int goto_dec_int = (int)floor(0.5+dec*(2147483648.0/M_PI));

      const long long int expected_duration = 10*1000000ULL;
      end_of_last_goto = now + expected_duration;

#ifdef PRINT_PERIODIC_POSITION
        std::cout << PrintRaInt(hour_angle_int)
                  << " FORBIDDEN, "
                     "expecting " << expected_duration << "us for "
                     "GOTO("
                  << "Ra:" << PrintRaInt(goto_ra_int)
                  << ",dec:" << PrintDecInt(goto_dec_int) << ')';
#endif
      gotoPosition(LX200RaDec(goto_ra_int,goto_dec_int));
    }
#ifdef PRINT_PERIODIC_POSITION
      std::cout << std:: endl;
#endif
  }
}


class TelescopeLX200::CommandGetPos : public TelescopeLX200::Command {
public:
  CommandGetPos(TelescopeLX200 &telescope) : Command(telescope) {}
  void execAsync(void) override {
    sendGetRaRqu();
  }
  void print(std::ostream &o) const override {o << "CommandGetPos";}
private:
  void sendGetRaRqu(void) {
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandGetPos::sendGetRaRqu: "
                 "sending :GR#" << std::endl;
    telescope.sendRqu(":GR#",4,std::bind(&TelescopeLX200::CommandGetPos::recvGetRaRsp,this));
  }
  int recvGetRaRsp(void) {
    if (telescope.recv_used < 8) return 0;
    const bool long_format = (telescope.recv_buf[5] == ':');
    if (long_format) {
      if (telescope.recv_used < 9) return 0;
    } else {
      if (telescope.recv_buf[5] != '.') {
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandGetPos::recvGetRaRsp: \""
                  << std::string(telescope.recv_buf,telescope.recv_used)
                  << "\": neigther long nor short format, ignoring response" << std::endl;
        return -1;
      }
    }
    if (telescope.recv_buf[7+long_format] != '#') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetRaRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": no '#' at end, ignoring response" << std::endl;
      return -1;
    }
    unsigned int hours;
    if (!ParseDecimal(telescope.recv_buf,2,hours)) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetRaRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad hours: no number, ignoring response" << std::endl;
      return -1;
    }
    if (hours >= 24) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetRaRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad hours: too big, ignoring response" << std::endl;
      return -1;
    }
    if (telescope.recv_buf[2] != ':') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetRaRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": no ':'l between hours and minutes, ignoring response" << std::endl;
      return -1;
    }
    unsigned int minutes;
    if (!ParseDecimal(telescope.recv_buf+3,2,minutes)) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetRaRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad minutes: no number, ignoring response" << std::endl;
      return -1;
    }
    if (minutes >= 60) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetRaRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad minutes: too big, ignoring response" << std::endl;
      return -1;
    }
    unsigned int seconds;
    if (long_format) {
      if (!ParseDecimal(telescope.recv_buf+6,2,seconds)) {
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandGetPos::recvGetRaRsp: \""
                  << std::string(telescope.recv_buf,telescope.recv_used)
                  << "\": bad seconds: no number, ignoring response" << std::endl;
        return -1;
      }
      if (seconds >= 60) {
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandGetPos::recvGetRaRsp: \""
                  << std::string(telescope.recv_buf,telescope.recv_used)
                  << "\": bad seconds: too big, ignoring response" << std::endl;
        return -1;
      }
    } else {
      if (!ParseDecimal(telescope.recv_buf+6,1,seconds)) {
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandGetPos::recvGetRaRsp: \""
                  << std::string(telescope.recv_buf,telescope.recv_used)
                  << "\": bad deciminutes: no number, ignoring response" << std::endl;
        return -1;
      }
      seconds += 6;
    }
    ra_dec.ra = (hours*60+minutes)*60+seconds;
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandGetPos::recvGetRaRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": ok: " << PrintRaSeconds(ra_dec.getRa()) << std::endl;
    if (long_format) {
      sendGetDecRqu();
      return 9;
    }
    sendPrecisionToggleRqu(std::bind(&TelescopeLX200::CommandGetPos::sendGetDecRqu,this));
    return 8;
  }    
  void sendPrecisionToggleRqu(std::function<void(void)> &&finished) {
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandGetPos::sendPrecisionToggleRqu: "
                 "sending :U#" << std::endl;
    telescope.sendMsg(":U#",3,std::move(finished));
  }
  void sendGetDecRqu(void) {
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandGetPos::sendGetDecRqu: "
                 "sending :GD#" << std::endl;
    telescope.sendRqu(":GD#",4,std::bind(&TelescopeLX200::CommandGetPos::recvGetDecRsp,this));
  }
  int recvGetDecRsp(void) {
    if (telescope.recv_used < 7) return 0;
      // new spec says '\'', but my scope sends ':'
    const bool long_format = (telescope.recv_buf[6] == ':' || telescope.recv_buf[6] == '\'');
    if (long_format) {
      if (telescope.recv_used < 10) return 0;
      if (telescope.recv_buf[9] != '#') {
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandGetPos::recvGetDecRsp: \""
                  << std::string(telescope.recv_buf,telescope.recv_used)
                  << "\": no '#' at end, ignoring response" << std::endl;
        return -1;
      }
    } else {
      if (telescope.recv_buf[6] != '#') {
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandGetPos::recvGetDecRsp: \""
                  << std::string(telescope.recv_buf,telescope.recv_used)
                  << "\": neigther long nor short format, ignoring response" << std::endl;
        return -1;
      }
    }
    ra_dec.dec_sign = (telescope.recv_buf[0] == '-');
    if (!ra_dec.dec_sign && (telescope.recv_buf[0] != '+')) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetDecRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": dec sign must be '+' or '-', ignoring response" << std::endl;
      return -1;
    }
    unsigned int degrees;
    if (!ParseDecimal(telescope.recv_buf+1,2,degrees)) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetDecRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad degrees: no number, ignoring response" << std::endl;
      return -1;
    }
      // spec: '*', scope: 223
    if (telescope.recv_buf[3] != (char)223 && telescope.recv_buf[3] != '*') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetDecRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": no degree symbol between degrees and minutes, ignoring response" << std::endl;
      return -1;
    }
    unsigned int minutes;
    if (!ParseDecimal(telescope.recv_buf+4,2,minutes)) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetDecRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad minutes: no number, ignoring response" << std::endl;
      return -1;
    }
    if (minutes >= 60) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetDecRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad minutes: too big, ignoring response" << std::endl;
      return -1;
    }
    unsigned int seconds;
    if (long_format) {
      if (!ParseDecimal(telescope.recv_buf+7,2,seconds)) {
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandGetPos::recvGetDecRsp: \""
                  << std::string(telescope.recv_buf,telescope.recv_used)
                  << "\": bad seconds: no number, ignoring response" << std::endl;
        return -1;
      }
      if (seconds >= 60) {
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandGetPos::recvGetDecRsp: \""
                  << std::string(telescope.recv_buf,telescope.recv_used)
                  << "\": bad seconds: too big, ignoring response" << std::endl;
        return -1;
      }
    } else {
      seconds = 0; // or rounding? 30?
    }
    ra_dec.dec = (degrees*60+minutes)*60+seconds;
    if (ra_dec.dec > 90*3600) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGetPos::recvGetDecRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad degrees: too big, ignoring response" << std::endl;
      return -1;
    }
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandGetPos::recvGetDecRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": ok: " << PrintDecSeconds(ra_dec.getDecSign(),ra_dec.getAbsDec()) << std::endl;
    telescope.positionReceived(ra_dec);
    telescope.get_pos_deadline.expires_from_now(boost::posix_time::microseconds(250000));
    telescope.get_pos_deadline.async_wait(
      [t = &telescope](const boost::system::error_code &e) {
        if (e) return; // timer was cancelled
        t->getPos();
      });
    if (long_format) {
      telescope.commandFinished();
      return 10;
    }
    sendPrecisionToggleRqu(std::bind(&TelescopeLX200::commandFinished,&telescope));
    return 7;
  }
private:
  LX200RaDec ra_dec;
};

void TelescopeLX200::getPos(void) {
  if (!next_command_get_pos) {
    next_command_get_pos = std::make_unique<CommandGetPos>(*this);
  }
  doSomething();
}

//------------------------------------------------------------------------




class TelescopeLX200::CommandGoto : public TelescopeLX200::Command {
public:
  CommandGoto(TelescopeLX200 &telescope) : Command(telescope) {}
  void set(const LX200RaDec &ra_dec) {
    CommandGoto::ra_dec = ra_dec;
  }
  void execAsync(void) override {
    sendDecRqu();
  }
  void print(std::ostream &o) const override {o << "CommandGoto" << ra_dec;}
private:
  void sendDecRqu(void) {
    buf[0] = ':';
    buf[1] = 'S';
    buf[2] = 'd';
    buf[3] = ra_dec.getDecSign() ? '-' : '+';
    unsigned int y = ra_dec.getAbsDec() / 10;
                             buf[11] = '0' + (char)(ra_dec.getAbsDec()-10*y);
    unsigned int x = y /  6; buf[10] = '0' + (char)(y- 6*x);
                             buf[ 9] = ':';
                 y = x / 10; buf[ 8] = '0' + (char)(x-10*y);
                 x = y /  6; buf[ 7] = '0' + (char)(y- 6*x);
                             buf[ 6] = '*';
                 y = x / 10; buf[ 5] = '0' + (char)(x-10*y);
                             buf[ 4] = '0' + (char)(     y);
                             buf[12] = '#';
                             buf[13] = '\0';
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandGoto" << ra_dec << "::sendDecRqu: "
                 "sending " << buf << std::endl;
    telescope.sendRqu(buf,13,std::bind(&TelescopeLX200::CommandGoto::recvDecRsp,this));
  }
  int recvDecRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGoto" << ra_dec << "::recvDecRsp \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                   "Sd: Bad Response" << std::endl;
      return -1;
    }
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandGoto" << ra_dec << "::recvDecRsp \""
              << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                 "ok" << std::endl;
    sendRaRqu();
    return 1;
  }
  void sendRaRqu(void) {
    buf[0] = ':';
    buf[1] = 'S';
    buf[2] = 'r';
    unsigned int y = ra_dec.getRa() / 10;
                             buf[10] = '0' + (char)(ra_dec.getRa()-10*y);
    unsigned int x = y /  6; buf[ 9] = '0' + (char)(y- 6*x);
                             buf[ 8] = ':';
                 y = x / 10; buf[ 7] = '0' + (char)(x-10*y);
                 x = y /  6; buf[ 6] = '0' + (char)(y- 6*x);
                             buf[ 5] = ':';
                 y = x / 10; buf[ 4] = '0' + (char)(x-10*y);
                             buf[ 3] = '0' + (char)(     y);
                             buf[11] = '#';
                             buf[12] = '\0';
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandGoto" << ra_dec << "::sendRaRqu: "
                 "sending " << buf << std::endl;
    telescope.sendRqu(buf,12,std::bind(&TelescopeLX200::CommandGoto::recvRaRsp,this));
  }
  int recvRaRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandGoto" << ra_dec << "::recvRaRsp \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                   "Sr: Bad Response" << std::endl;
      return -1;
    }
    sendGotoRqu();
    return 1;
  }
  void sendGotoRqu(void) {
    buf[0] = ':';
    buf[1] = 'M';
    buf[2] = 'S';
    buf[3] = '#';
    std::cout << PrintTime() << " "
                 "TelescopeLX200::CommandGoto" << ra_dec << "::sendGotoRqu: "
                 "sending :MS#" << std::endl;
    telescope.sendRqu(buf,4,std::bind(&TelescopeLX200::CommandGoto::recvGotoRsp,this));
  }
  int recvGotoRsp(void) {
    switch (telescope.recv_buf[0]) {
      case '0':
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandGoto" << ra_dec << "::recvGotoRsp \""
                  << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                     "GOTO ok."
                  << std::endl;
        telescope.commandFinished();
        return 1;
      case '1':
      case '2':
        for (int i=1;i<telescope.recv_used;i++) {
          if (telescope.recv_buf[i] == '#') {
            std::cout << PrintTime() << " "
                         "TelescopeLX200::CommandGoto" << ra_dec << "::recvGotoRsp: "
                         "below " << ((telescope.recv_buf[0]=='1')?"Horizon":"Higher Limit")
                      << ", \""
                      << std::string(telescope.recv_buf,telescope.recv_used) << '"' << std::endl;
            telescope.commandFinished();
            return (i+1);
          }
        }
        return 0;
      default:
        std::cout << PrintTime() << " "
                     "TelescopeLX200::CommandGoto" << ra_dec << "::recvGotoRsp("
                  << std::string(telescope.recv_buf,telescope.recv_used) << "): "
                     "MS: Bad Response" << std::endl;
        return -1;
    }
  }
private:
  LX200RaDec ra_dec;
  char buf[14];
};

void TelescopeLX200::gotoPosition(const LX200RaDec &ra_dec) {
  std::cout << PrintTime() << " "
               "TelescopeLX200::gotoPosition" << ra_dec << std::endl;
  if (!next_command_goto) {
    next_command_goto = std::make_unique<CommandGoto>(*this);
  }
  next_command_goto->set(ra_dec);
  doSomething();
}

void TelescopeLX200::gotoPosition(const unsigned int ra_int_j2000,const int dec_int_j2000) {
//    if (dec_int_j2000 < -0x40000000 || dec_int_j2000 > 0x40000000) abort();

  const Vector<double,3> v0 = PolarToRect( ra_int_j2000 * (M_PI/2147483648.0),
                                          dec_int_j2000 * (M_PI/2147483648.0) );
  const Vector<double,3> v = v0*precession_matrix;
  double ra,dec;
  RectToPolar(v,ra,dec);
  if (dec < -0.5*M_PI || dec > 0.5*M_PI) abort();
  LX200RaDec ra_dec( (unsigned int)floor(0.5+ ra*(2147483648.0/M_PI)),
                                (int)floor(0.5+dec*(2147483648.0/M_PI)) );
  gotoPosition(ra_dec);
}


//------------------------------------------------------------------------





class TelescopeLX200::CommandMove : public TelescopeLX200::Command {
public:
  CommandMove(TelescopeLX200 &telescope) : Command(telescope) {}
  void accumulate(short int horz,short int vert,
                  unsigned int micros) {
    if (horz != -0x8000) CommandMove::horz = horz;
    if (vert != -0x8000) CommandMove::vert = vert;
    CommandMove::micros = micros;
    std::cout << "CommandMove::accumulate: "
              << CommandMove::horz << '/' << CommandMove::vert << ": " << CommandMove::micros
              << std::endl;
  }
  static short int Rescale9(short int x) {
      // range -9..+9
//      return (x < 0) ? -(-x*9+0x4000) >> 15) : (x*9+0x4000) >> 15);
      // 2*32767/19=3449.16
    return (x < 0) ? -((-x + 1725) / 3450) : ((x + 1725) / 3450);
  }
  void execAsync(void) override {
    if (horz == -0x8000) {
      horz = telescope.curr_horz;
    } else {
        // range -9..+9
      horz = Rescale9(horz);
    }
    if (vert == -0x8000) {
      vert = telescope.curr_vert;
    } else {
      vert = Rescale9(vert);
    }
    speed = std::max(std::abs(horz),std::abs(vert));
    if (std::abs(horz) < speed) horz = 0;
    if (std::abs(vert) < speed) vert = 0;
    std::cout << "CommandMove::execAsync: "
              << horz << '/' << vert
              << std::endl;
    if (speed != 0) {
      telescope.move_deadline.expires_from_now(boost::posix_time::microseconds(micros));
      telescope.move_deadline.async_wait(
        [t = &telescope](const boost::system::error_code &e) {
          std::cout << "CommandMove::move_deadline_l: " << e.message() << std::endl;
          if (e) return; // timer was cancelled
          std::cout << "CommandMove::move_deadline_l: stopping movement" << std::endl;
          t->move(0,0,0);
        });
    }
    if (telescope.curr_horz == horz && telescope.curr_vert == vert) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandMove::execAsync: already stopped, nothing to do." << std::endl;
        // nothing to do
      telescope.commandFinished();
      return;
    }
    sendSpeedRqu();
  }
  void print(std::ostream &o) const override {
    o << "CommandMove(" << horz << ',' << vert << ';' << micros << ')';
  }
private:
  void sendSpeedRqu(void) {
    if (speed == 0) {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandMove::sendSpeedRqu: stop" << std::endl;
      telescope.move_deadline.cancel();
      sendStopRqu();
      return;
    }
    if (telescope.curr_speed == speed) {
      sendHorzRqu();
      return;
    }
    if (telescope.curr_speed == 0) {
        // start moving, TODO: query movement speed for later resetting
    }
    buf[0] = ':';
    buf[1] = 'S';
    buf[2] = 'R';
    buf[3] = '0'+(char)speed;
    buf[4] = '#';
    telescope.sendRqu(buf,5,std::bind(&TelescopeLX200::CommandMove::recvSpeedRsp,this));
  }
  int recvSpeedRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandMove(" << horz << ',' << vert << ")::recvSpeedRsp("
                << std::string(telescope.recv_buf,telescope.recv_used) << "): "
                   "Bad Response" << std::endl;
      return -1;
    }
    telescope.curr_speed = speed;
    sendHorzRqu();
    return 1;
  }
  void sendStopRqu(void) {
    telescope.sendRqu(":q#",3,std::bind(&TelescopeLX200::CommandMove::recvStopRsp,this));
  }
  int recvStopRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandMove(" << horz << ',' << vert << ")::recvStopRsp("
                << std::string(telescope.recv_buf,telescope.recv_used) << "): "
                   "Bad Response" << std::endl;
      return -1;
    }
    telescope.curr_speed = 0;
    telescope.curr_horz = 0;
    telescope.curr_vert = 0;
    telescope.commandFinished();
    return 1;
  }
  void sendHorzRqu(void) {
    if (telescope.curr_horz == horz) {
      sendVertRqu();
      return;
    }
    buf[0] = ':';
    buf[3] = '#';
    if (horz == 0) {
      buf[1] = 'q';
      buf[2] = 'R';
      telescope.sendRqu(buf,4,std::bind(&TelescopeLX200::CommandMove::recvHorzStopRsp,this));
    } else {
      buf[1] = 'm';
        // horz > 0: right, east
      buf[2] = (horz < 0) ? 'w' : 'e';
      telescope.sendMsg(buf,4,std::bind(&TelescopeLX200::CommandMove::horzRquFinished,this));
    }
  }
  void horzRquFinished(void) {
    telescope.curr_horz = horz;
    sendVertRqu();
  }
  int recvHorzStopRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandMove(" << horz << ',' << vert << ")::recvHorzStopRsp("
                << std::string(telescope.recv_buf,telescope.recv_used) << "): "
                   "Bad Response" << std::endl;
      return -1;
    }
    horzRquFinished();
    return 1;
  }
  void sendVertRqu(void) {
    if (telescope.curr_vert == vert) {
      telescope.commandFinished();
      return;
    }
    buf[0] = ':';
    buf[3] = '#';
    if (vert == 0) {
      buf[1] = 'q';
      buf[2] = 'D';
      telescope.sendRqu(buf,4,std::bind(&TelescopeLX200::CommandMove::recvVertStopRsp,this));
    } else {
      buf[1] = 'm';
        // vert > 0: down, south
      buf[2] = (vert < 0) ? 's' : 'n'; // :mn# actually moves south
      buf[5] = '\0';
std::cout << "vertRqu: " << buf << std::endl;
      telescope.sendMsg(buf,4,std::bind(&TelescopeLX200::CommandMove::vertRquFinished,this));
    }
  }
  void vertRquFinished(void) {
    telescope.curr_vert = vert;
    telescope.commandFinished();
  }
  int recvVertStopRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeLX200::CommandMove(" << vert << ',' << vert << ")::recvVertStopRsp("
                << std::string(telescope.recv_buf,telescope.recv_used) << "): "
                   "Bad Response" << std::endl;
      return -1;
    }
    vertRquFinished();
    return 1;
  }
private:
  unsigned int micros = 0;
  short int horz = 0;
  short int vert = 0;
  short int speed = 0;
  char buf[5];
};

void TelescopeLX200::move(short int horz,short int vert,
          unsigned int validity_micros) {
  if (!next_command_move) {
    next_command_move = std::make_unique<CommandMove>(*this);
  }
  next_command_move->accumulate(horz,vert,validity_micros);
  next_command_goto.reset();
  doSomething();
}


//------------------------------------------------------------------------




class TelescopeLX200::CommandGuide  : public TelescopeLX200::Command {
public:
  CommandGuide(TelescopeLX200 &telescope) : Command(telescope),d_ra(0),d_dec(0) {}
  void accumulate(int d_ra_micros,int d_dec_micros) {
    if (d_ra_micros != (-0x7FFFFFFF-1)) d_ra += d_ra_micros;
    if (d_dec_micros != (-0x7FFFFFFF-1)) d_dec = d_dec_micros;
  }
  void execAsync(void) override {
      // convert micros->millis
    d_ra /= 1000;
    d_dec /= 1000;
    sendRaRqu();
  }
  void print(std::ostream &o) const override {
    o << "CommandGuide(" << d_ra << ',' << d_dec << ')';
  }
private:
  void sendRaRqu(void) {
    if (d_ra == 0) {
      sendDecRqu();
      return;
    }
    buf[0] = ':';
    buf[1] = 'M';
    buf[2] = (d_ra < 0) ? 'w' : 'e';
    FillDecimal(buf+3,5,std::abs(d_ra));
    buf[8] = '#';
    telescope.sendMsg(buf,9,std::bind(&TelescopeLX200::CommandGuide::raRquFinished,this));
  }
  void raRquFinished(void) {
    d_ra = 0;
    sendDecRqu();
  }
  void sendDecRqu(void) {
    if (d_dec == 0) {
      telescope.commandFinished();
      return;
    }
    buf[0] = ':';
    buf[1] = 'M';
    buf[2] = (d_dec < 0) ? 'n' : 's';
    FillDecimal(buf+3,5,std::abs(d_dec));
    buf[8] = '#';
    telescope.sendMsg(buf,9,std::bind(&TelescopeLX200::CommandGuide::decRquFinished,this));
  }
  void decRquFinished(void) {
    d_dec = 0;
    telescope.commandFinished();
  }
private:
  int d_ra;
  int d_dec;
  char buf[9];
};


void TelescopeLX200::guide(int d_ra_micros,int d_dec_micros) {
  if (!next_command_guide) {
    next_command_guide = std::make_unique<CommandGuide>(*this);
  }
  next_command_guide->accumulate(d_ra_micros,d_dec_micros);
  doSomething();
}


//------------------------------------------------------------------------





void TelescopeLX200::init(unsigned int drain_micros) {
  get_pos_deadline.cancel();
  if (!next_command_init) {
    next_command_init = std::make_unique<CommandInit>(*this,drain_micros);
  }
  next_command_get_loc.reset();
  next_command_get_pos.reset();
  next_command_goto.reset();
  next_command_move.reset();
  next_command_guide.reset();
  commandFinished();
}

void TelescopeLX200::commandFinished(void) {
//    std::cout << "commandFinished" << std::endl;
  command_deadline.cancel();
  curr_command.reset();
  doSomething();
}

void TelescopeLX200::doSomething(void) {
//    std::cout << "doSomething start" << std::endl;
  if (curr_command) return;
  if (next_command_init) {
//      std::cout << "doSomething: init" << std::endl;
    curr_command = std::move(next_command_init);
  } else if (next_command_get_loc) {
//      std::cout << "doSomething: get_loc" << std::endl;
    curr_command = std::move(next_command_get_loc);
  } else if (next_command_get_pos) {
//      std::cout << "doSomething: get_pos" << std::endl;
    curr_command = std::move(next_command_get_pos);
  } else if (next_command_move) {
//      std::cout << "doSomething: move" << std::endl;
    curr_command = std::move(next_command_move);
  } else if (next_command_goto) {
//      std::cout << "doSomething: goto" << std::endl;
    curr_command = std::move(next_command_goto);
  } else if (next_command_guide) {
//      std::cout << "doSomething: guide" << std::endl;
    curr_command = std::move(next_command_guide);
  } else {
    // nothing to do
//      std::cout << "doSomething: nothing" << std::endl;
  }
  if (curr_command) {
//      std::cout << "doSomething: sceduling curr_command" << std::endl;
    curr_command->execAsync();
  }
    // execAsync() might have called call commandFinished()
    // which resets curr_command and calls doSomething(): harmless.
    // But then curr_command might be invalid, so check again:
  if (curr_command) {
    command_deadline.expires_from_now(boost::posix_time::microseconds(curr_command->getTimeoutMicros()));
    command_deadline.async_wait(
      [this](const boost::system::error_code &error) {
        if (error) {
//          std::cout << PrintTime() << " "
//                       "TelescopeLX200::doSomething::timeout_l: " << error.message() << std::endl;
          return; // probably cancelled
        }
        std::cout << PrintTime() << " "
                     "TelescopeLX200::doSomething::timeout_l: " << *curr_command << " timeout" << std::endl;
        serial.cancel();
        doSomething();
      });
  }
//    std::cout << "doSomething end" << std::endl;
}



