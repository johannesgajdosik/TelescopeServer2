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







struct IOptronRaDec {
  IOptronRaDec(void) {}
  IOptronRaDec(unsigned int ra,bool dec_sign,unsigned int dec)
    : ra(ra),dec_sign(dec_sign),dec(dec) {
    if (dec > 90*3600*100) abort();
  }
  unsigned int getRa(void) const {return ra;}
  bool getDecSign(void) const {return dec_sign;}
  unsigned int getAbsDec(void) const {return dec;}
  int getDec(void) const {return (dec_sign ? -dec : dec);}
  IOptronRaDec(unsigned int ra_int,int dec_int) {
    set(ra_int,dec_int);
  }
  void set(unsigned int ra_int,int dec_int) {
    ra  = ((unsigned int)((                                                      ra_int  * (24*3600*1000ULL) + 0x80000000ULL) >> 32));
    dec_sign = (dec_int < 0);
    dec = ((unsigned int)((  (dec_sign ? -(unsigned int)dec_int : (unsigned int)dec_int) * (360*3600*100ULL) + 0x80000000ULL) >> 32));
// debugging:
//    unsigned int abs_dec_int = (dec_sign ? -(unsigned int)dec_int : (unsigned int)dec_int);
//    unsigned long long int d2 = abs_dec_int * (360*3600*100ULL) + 0x80000000ULL;
//    unsigned int d2_32 = d2 >> 32;
//    if (d2_32 != dec) abort();
//    if (dec > 90*3600*100) abort();
  }
  unsigned int getRaInt(void) const {
    return ((((unsigned long long int)ra)<<32)+(12*3600*1000ULL))/(24*3600*1000ULL);
  }
  double getRaRad(void) const {return ra*(M_PI/(12*3600*1000.0));}
  int getDecInt(void) const {
    int dec_int = ((((unsigned long long int)dec)<<32)+(180*3600*100ULL))/(360*3600*100ULL);
    if (dec_sign) dec_int = -dec_int;
    return dec_int;
  }
  double getDecRad(void) const {
    const double dec_rad = dec*(M_PI/(180*3600*100.0));
    return (dec_sign ? -dec_rad : dec_rad);
  }
  unsigned int ra;
  bool dec_sign;
  unsigned int dec;
};

static
std::ostream &operator<<(std::ostream &o,const IOptronRaDec &x) {
  o << "(Ra:" << PrintRaMilliseconds(x.getRa())
    << ",Dec:" << PrintDecCentiseconds(x.getDecSign(),x.getAbsDec()) << ')';
  return o;
}









class TelescopeIOptron : public Telescope {
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
    Command(TelescopeIOptron &telescope) : telescope(telescope) {}
    TelescopeIOptron &telescope;
  };
  TelescopeIOptron(const TelescopeIOptron&) = delete;
  TelescopeIOptron &operator=(const TelescopeIOptron&) = delete;
private:
  TelescopeIOptron(const std::string &args,
                   OpenedClosedFunction &&opened_closed,
                   PositionFunction &&announce_position,
                   boost::asio::io_context &io_context);
  ~TelescopeIOptron(void);
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
  void locationReceived(unsigned int longitude_ioptron,unsigned int latitude_p90_ioptron);
  class CommandGetLoc;
  void getLoc(void);

    // telescope pointing position in the sky
  void positionReceived(const IOptronRaDec &ra_dec);
  class CommandGetPos;
  void getPos(void);

    // Goto
  class CommandGoto;
  void gotoPosition(const IOptronRaDec &ra_dec);
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

  char recv_buf[32];
  int recv_used;
  short int curr_horz = 0;
  short int curr_vert = 0;
  short int curr_speed = 0;

  bool telescope_online = false;
  bool language_v3 = false;
};

static inline
std::ostream &operator<<(std::ostream &o,const TelescopeIOptron::Command &c) {
  c.print(o);
  return o;
}

static const bool telescope_IOptron_registered
  = Telescope::RegisterCreationFunction("iOptron",TelescopeIOptron::Create);






















//-----------------balcony railing---------------------------------------------
/*
before crash:

Dec: +79d12m21s03, Az: 12h25m29s404, Alt: +38d27m40s11, h: 10h09m24s091
Dec: +73d00m03s80, Az: 12h30m55s421, Alt: +32d15m19s21, h: 10h28m24s026
Dec: +70d18m54s49, Az: 12h31m50s393, Alt: +29d30m30s88, h: 10h36m08s624
Dec: +65d19m08s37, Az: 12h31m04s374, Alt: +24d18m54s04, h: 10h51m22s386
Dec: +59d41m39s85, Az: 12h28m16s972, Alt: +18d27m52s71, h: 11h06m29s195
Dec: +55d40m33s29, Az: 12h23m47s203, Alt: +14d15m02s26, h: 11h18m58s085
Dec: +48d37m37s67, Az: 12h18m39s837, Alt: +07d02m14s87, h: 11h31m56s149


Dec: +84d59m42s13, Az: 12h14m38s792, Alt: +43d54m22s41, h: 09h52m40s614
Dec: +80d06m42s54, Az: 12h24m45s401, Alt: +39d23m21s27, h: 10h03m53s769
Dec: +75d02m41s96, Az: 12h30m41s344, Alt: +34d23m59s31, h: 10h18m55s419
Dec: +70d01m02s76, Az: 12h33m05s248, Alt: +29d16m37s69, h: 10h33m49s316
Dec: +65d03m04s10, Az: 12h32m42s093, Alt: +24d07m27s23, h: 10h48m19s270
Dec: +59d55m01s64, Az: 12h30m37s680, Alt: +18d47m10s97, h: 11h01m41s846
Dec: +55d04m14s65, Az: 12h27m26s927, Alt: +13d45m26s56, h: 11h13m13s342
Dec: +50d01m15s00, Az: 12h21m22s188, Alt: +08d29m48s41, h: 11h27m02s355


Dec: +85d00m01s53, Az: 12h14m46s890, Alt: +43d55m33s88, h: 09h51m15s513
Dec: +89d07m16s23, Az: 12h02m58s365, Alt: +47d29m46s81, h: 09h40m37s007

Dec: +45d00m16s18, Az: 12h16m21s044, Alt: +03d21m42s37, h: 11h36m53s700
Dec: +40d00m55s44, Az: 12h12m01s402, Alt: -01d41m45s66, h: 11h44m18s171
Dec: +35d01m58s48, Az: 12h12m54s422, Alt: -06d40m25s11, h: 11h44m20s398
Dec: +29d59m30s48, Az: 12h27m29s408, Alt: -11d28m38s44, h: 11h28m52s401
Dec: +25d06m47s99, Az: 12h42m50s033, Alt: -15d57m31s29, h: 11h14m28s989
Dec: +19d59m53s49, Az: 12h59m03s386, Alt: -20d31m13s73, h: 11h01m08s785
Dec: +15d00m29s29, Az: 13h15m58s085, Alt: -24d48m34s45, h: 10h48m46s005
Dec: +09d59m59s54, Az: 13h32m37s868, Alt: -29d02m12s21, h: 10h38m15s684
*/

struct DecHA {
  int dec_int;
  int ha_int;
};

static constexpr int dec_ha_limits_size = 13;
static constexpr DecHA dec_ha_limits[dec_ha_limits_size] = {
  { (int)(30*0x80000000LL/180), (int)((11*60+28)*0x80000000LL/(12*60)) },
  { (int)(35*0x80000000LL/180), (int)((11*60+44)*0x80000000LL/(12*60)) },
  { (int)(40*0x80000000LL/180), (int)((11*60+44)*0x80000000LL/(12*60)) },
  { (int)(45*0x80000000LL/180), (int)((11*60+37)*0x80000000LL/(12*60)) },
  { (int)(50*0x80000000LL/180), (int)((11*60+27)*0x80000000LL/(12*60)) },
  { (int)(55*0x80000000LL/180), (int)((11*60+13)*0x80000000LL/(12*60)) },
  { (int)(60*0x80000000LL/180), (int)((11*60+ 2)*0x80000000LL/(12*60)) },
  { (int)(65*0x80000000LL/180), (int)((10*60+48)*0x80000000LL/(12*60)) },
  { (int)(70*0x80000000LL/180), (int)((10*60+34)*0x80000000LL/(12*60)) },
  { (int)(75*0x80000000LL/180), (int)((10*60+19)*0x80000000LL/(12*60)) },
  { (int)(80*0x80000000LL/180), (int)((10*60+04)*0x80000000LL/(12*60)) },
  { (int)(85*0x80000000LL/180), (int)(( 9*60+53)*0x80000000LL/(12*60)) },
  { (int)(90*0x80000000LL/180), (int)(( 9*60+41)*0x80000000LL/(12*60)) }
};

static constexpr
bool ForbiddenDecHA(int dec_int,int ha_int) {
  if (dec_int <= dec_ha_limits[0].dec_int) return false;
  for (int i=1;i<dec_ha_limits_size;i++) {
    if (dec_int <= dec_ha_limits[i].dec_int) {
         // dec_ha_limits[i-1].dec_int <= dec_int < dec_ha_limits[i].dec_int
      return ( ha_int * (long long int)(dec_ha_limits[i].dec_int - dec_ha_limits[i-1].dec_int)
             > dec_ha_limits[i-1].ha_int * (long long int)(dec_ha_limits[i].dec_int - dec_int)
             + dec_ha_limits[i].ha_int * (long long int)(dec_int - dec_ha_limits[i-1].dec_int) );
    }
  }
  return true;
}
//------------------------------------------------------------------------------





//-----construction/destruction---------------------------------------------------


Telescope::Ptr
TelescopeIOptron::Create(const std::string &args,
                         OpenedClosedFunction &&opened_closed,
                         PositionFunction &&announce_position,
                         boost::asio::io_context &io_context) {
  return new TelescopeIOptron(args,
                              std::move(opened_closed),
                              std::move(announce_position),
                              io_context);
}

TelescopeIOptron::TelescopeIOptron(const std::string &args,
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
  io_context.dispatch(std::bind(&TelescopeIOptron::initialize,this));
  std::cout << PrintTime() << " "
               "TelescopeIOptron::TelescopeIOptron(" << args << ')' << std::endl;
}

TelescopeIOptron::~TelescopeIOptron(void) {
  if (telescope_online) {
    opened_closed(false);
    telescope_online = false;
  }
  std::cout << PrintTime() << " "
               "TelescopeIOptron::~TelescopeIOptron" << std::endl;
}

//------------------------------------------------------------------------------




//-------------------RS232 communication------------------------------------


void TelescopeIOptron::recvRsp(std::function<int(void)> &&rsp_data_received) {
  const int to_read = (int)(sizeof(recv_buf)-1) - recv_used;
  if (to_read <= 0) {
    std::cout << PrintTime() << " "
                 "TelescopeIOptron::recvRsp: recv_buf full, discarding response"
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
                       "TelescopeIOptron::recvRsp::l: "
                       "IO ERROR" << std::endl;
          initialize();
          return;
        }
        if (error == boost::asio::error::eof) {
          std::cout << PrintTime() << " "
                       "TelescopeIOptron::recvRsp::l: "
                       "EOF" << std::endl;
          initialize();
          return;
        }
        if (error == boost::asio::error::operation_aborted) {
//          std::cout << PrintTime() << " "
//                       "TelescopeIOptron::recvRsp::l: "
//                       "canceled" << std::endl;
          init();
          return;
        }
        std::cout << PrintTime() << " "
                     "TelescopeIOptron::recvRsp::l: read failed: "
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

void TelescopeIOptron::sendRqu(const char *data,int size,
                               std::function<int(void)>  &&rsp_data_received) {
  sendMsg(data,size,
          [this,f = std::move(rsp_data_received)]() mutable {
            recv_used = 0;
            std::function<int(void)> f2;
            f2.swap(f);
            recvRsp(std::move(f2));
          });
}

void TelescopeIOptron::sendMsg(const char *data,int size,
                               std::function<void(void)> &&finished) {
  serial.async_write_some(
    boost::asio::buffer(data,size),
    [this,data,size,f = std::move(finished)] (
        const boost::system::error_code &error,std::size_t written) mutable {
      if (error) {
        if (error == boost::system::errc::io_error) {
          std::cout << PrintTime() << " "
                       "TelescopeIOptron::sendMsg::l: "
                       "IO ERROR" << std::endl;
          initialize();
          return;
        }
        std::cout << PrintTime() << " "
                     "TelescopeIOptron::sendMsg::l: write failed: "
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

void TelescopeIOptron::initPrecessionMatrix(void) {
  const double jd_minus_j2000 = TimeToJD_minus_J2000(GetNow());
  const double t = jd_minus_j2000 * (1.0 / 365250.0);
  Lea406aPrecMat(t,precession_matrix);
}

void TelescopeIOptron::initialize(void) {
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
                 "TelescopeIOptron::initialize: serial port::close failed: "
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
                 "TelescopeIOptron::initialize: serial port::open(" << args << ") failed: "
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
               "TelescopeIOptron::initialize: serial port " << args << " opened" << std::endl;
  initPrecessionMatrix();
  init(0);
}


static const char *mount_info[] = {
  "0010", "160610", "Cube II or Cube Pro EQ mode",
  "0011", "161028", "SmartEQ Pro+",
  "0025", "170518", "CEM25(/P)",
  "0026", "170518", "CEM25-EC",
  "0030", "161101", "iEQ30 Pro",
  "0040", "181018", "CEM40",
  "0041", "181018", "CEM40-EC",
  "0045", "161101", "iEQ45 Pro EQ mode",
  "0046", "161101", "iEQ45 Pro AA mode", 
  "0060", "161101", "CEM60",
  "0061", "161101", "CEM60-EC",
  "0120", "171001", "CEM120",
  "0121", "171001", "CEM120-EC",
  "0122", "171001", "CEM120-EC2",
  "5010", "160610", "Cube II or Cube Pro AA mode",
  "5035", "170410", "AZ Mount Pro",
  0
};

static
bool GetMountInfoString(const char rsp[4],const char *&min_firmware,const char *&name) {
  const char **p = mount_info;
  while (*p) {
    if (0 == memcmp(p[0],rsp,4)) {
      min_firmware = p[1];
      name = p[2];
      return true;
    }
    p += 3;
  }
  return false;
}

class TelescopeIOptron::CommandInit : public TelescopeIOptron::Command {
public:
  CommandInit(TelescopeIOptron &telescope,unsigned int drain_micros)
    : Command(telescope),drain_deadline(telescope.serial.get_executor()),
      drain_micros(drain_micros) {}
  void execAsync(void) override {
    sendStopRqu();
  }
  void print(std::ostream &o) const override {
    o << "CommandInit(" << drain_micros << ')';
  }
private:
  unsigned int getTimeoutMicros(void) const override {return 500000 + drain_micros;}
  void sendStopRqu(void) {
    if (telescope.telescope_online) {
      telescope.telescope_online = false;
      telescope.opened_closed(false);
    }
    if (drain_micros == 0) {
      sendMountInfoRqu();
      return;
    }
    drain_deadline.expires_from_now(boost::posix_time::microseconds(drain_micros));
    drain_deadline.async_wait(
      [t = &telescope](const boost::system::error_code &error) {
//        std::cout << PrintTime() << " "
//                     "TelescopeIOptron::CommandInit::sendStopRqu::l: "
//                     "draining deadline reached, " << error.message()
//                  << std:: endl;
          // cancel serial read/write requests regardless of error:
        boost::system::error_code ec;
        t->serial.cancel(ec);
        if (ec) {
          std::cout << PrintTime() << " "
                       "TelescopeIOptron::CommandInit::sendStopRqu::l: "
                       "serial.cancel failed: " << ec.message()
                    << std:: endl;
        }
      });
    telescope.sendMsg(":q#",3,std::bind(&TelescopeIOptron::CommandInit::recvStopRsp,this));
  }
  void recvStopRsp(void) {
    telescope.serial.async_read_some(
      boost::asio::buffer(telescope.recv_buf,sizeof(telescope.recv_buf)-1),
      [this](const boost::system::error_code &error,std::size_t size) {
        if (error) {
          if (error == boost::system::errc::io_error) {
            std::cout << PrintTime() << " "
                         "TelescopeIOptron::CommandInit::recvStopRsp::l: "
                         "IO ERROR" << std::endl;
            telescope.initialize();
            return;
          }
          if (error == boost::asio::error::eof) {
            std::cout << PrintTime() << " "
                         "TelescopeIOptron::CommandInit::recvStopRsp::l: "
                         "EOF" << std::endl;
            telescope.initialize();
            return;
          }
          if (error == boost::asio::error::operation_aborted) {
//            std::cout << PrintTime() << " "
//                         "TelescopeIOptron::CommandInit::recvStopRsp::l: "
//                         "draining finished"
//                      << std:: endl;
            telescope.recv_used = 0;
            sendMountInfoRqu();
            return;
          }
          std::cout << PrintTime() << " "
                       "TelescopeIOptron::CommandInit::recvStopRsp::l: "
                       "unexpected error: " << error.message() << std::endl;
          telescope.initialize();
          return;
          
        }
        telescope.recv_buf[size] = '\0';
        std::cout << PrintTime() << " "
                     "TelescopeIOptron::CommandInit::recvStopRsp::l: "
                     "received \"" << telescope.recv_buf << "\", draining..."
                  << std::endl;
        recvStopRsp();
      });
  }
  void sendMountInfoRqu(void) {
    telescope.sendRqu(":MountInfo#",11,std::bind(&TelescopeIOptron::CommandInit::recvMountInfoRsp,this));
  }
  int recvMountInfoRsp(void) {
    if (telescope.recv_used < 4) return 0;
    const char *mount_name;
    if (!GetMountInfoString(telescope.recv_buf,min_firmware,mount_name)) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandInit::recvMountInfoRsp: Unknown Mount: \""
                << std::string(telescope.recv_buf,4) << '"' << std::endl;
      return -1;
    }
    std::cout << PrintTime() << " "
                 "TelescopeIOptron::CommandInit::recvMountInfoRsp: " << mount_name << " detected" << std::endl;
    if (!strncmp(telescope.recv_buf,"0120",4) || !strncmp(telescope.recv_buf,"0121",4) ||
        !strncmp(telescope.recv_buf,"0122",4)) {
        // TODO: also support CEM70 with language_v3
      telescope.language_v3 = true;
    }
    sendFirmware1Rqu();
    return 1;
  }
  void sendFirmware1Rqu(void) {
    telescope.sendRqu(":FW1#",5,std::bind(&TelescopeIOptron::CommandInit::recvFirmware1Rsp,this));
  }
  int recvFirmware1Rsp(void) {
    if (telescope.recv_used < 13) return 0;
    if (telescope.recv_buf[12] != '#') {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandInit::recvMountInfoRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used) << ": "
                   "Bad Response" << std::endl;
      return -1;
    }
    std::cout << PrintTime() << " "
                 "TelescopeIOptron::CommandInit::recvFirmware1Rsp: mainboard FW: "
              << std::string(telescope.recv_buf,6)
              << ", hand controller FW: "
              << std::string(telescope.recv_buf+6,6)
              << std::endl;
    if (memcmp(telescope.recv_buf,min_firmware,6) < 0) {
      std::cout << "TelescopeIOptron::CommandInit::recvFirmware1Rsp: "
                   "Mainboard firmware is too old, minimum is 20" << min_firmware << std::endl;
      return -1;
    }
    sendFirmware2Rqu();
    return 1;
  }
  void sendFirmware2Rqu(void) {
    telescope.sendRqu(":FW2#",5,std::bind(&TelescopeIOptron::CommandInit::recvFirmware2Rsp,this));
  }
  int recvFirmware2Rsp(void) {
    if (telescope.recv_used < 13) return 0;
    if (telescope.recv_buf[12] != '#') {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandInit::recvMountInfoRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used) << ": "
                   "Bad Response" << std::endl;
      return -1;
    }
    std::cout << PrintTime() << " "
                 "TelescopeIOptron::CommandInit::recvFirmware2Rsp: RA motor board FW: "
              << std::string(telescope.recv_buf,6)
              << ", Dec motor board FW: "
              << std::string(telescope.recv_buf+6,6)
              << std::endl;
    telescope.getLoc();
    if (!telescope.telescope_online) {
      telescope.telescope_online = true;
      telescope.opened_closed(true);
    }
    telescope.commandFinished();
    return 1;
  }
private:
  boost::asio::deadline_timer drain_deadline;
  const unsigned int drain_micros;
  const char *min_firmware;
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

void TelescopeIOptron::locationReceived(unsigned int longitude_ioptron,unsigned int latitude_p90_ioptron) {
  const unsigned int new_longitude_int =
       // revert longitude sign: east is positive for IOptron but negative for my code
    -(unsigned int)(((((unsigned long long int)longitude_ioptron) << 32) + 180*3600ULL) / (360*3600ULL));
  const int new_latitude_int =
    ((int)(((((unsigned long long int)latitude_p90_ioptron) << 32) + 180*3600ULL) / (360*3600ULL))) - 0x40000000;

  if (longitude_int != new_longitude_int ||
      latitude_int != new_latitude_int) {
    longitude_int = new_longitude_int;
    latitude_int = new_latitude_int;
    std::cout << PrintTime() << " "
                 "TelescopeIOptron::locationReceived: new location ("
              << PrintRaInt(longitude_int) << ','
              << PrintDecInt(latitude_int) << ')' << std::endl;
    SetLongLat(longitude_int*(M_PI/2147483648.0),
               latitude_int*(M_PI/2147483648.0),
               geographic_pos_orientation);
  }
}

class TelescopeIOptron::CommandGetLoc : public TelescopeIOptron::Command {
public:
  CommandGetLoc(TelescopeIOptron &telescope) : Command(telescope) {}
  void execAsync(void) override {
    sendGetLocRqu();
  }
  void print(std::ostream &o) const override {o << "CommandGetLoc";}
private:
  void sendGetLocRqu(void) {
    telescope.sendRqu(":GLS#",5,std::bind(&TelescopeIOptron::CommandGetLoc::recvGetLocRsp,this));
  }
  int recvGetLocRsp(void) {
    if (telescope.recv_used < 20) return 0;
    if (telescope.recv_buf[19] != '#') {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetLoc::recvGetLocRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": no '#' at end, ignoring response" << std::endl;
      return -1;
    }
    const bool longitude_sign = (telescope.recv_buf[0] == '-');
    if (!longitude_sign && (telescope.recv_buf[0] != '+')) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetLoc::recvGetLocRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": longitude sign must be '+' or '-', ignoring response" << std::endl;
      return -1;
    }
    unsigned int abs_longitude;
    if (!ParseDecimal(telescope.recv_buf+1,6,abs_longitude)) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetLoc::recvGetLocRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad longitude: no number, ignoring response" << std::endl;
      return -1;
    }
    if (abs_longitude > 180*3600) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetLoc::recvGetLocRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad longitude: abs(longitude)>" << (180*3600) << ", ignoring response"  << std::endl;
      return -1;
    }
    unsigned int latitude_p90;
    if (!ParseDecimal(telescope.recv_buf+7,6,latitude_p90)) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetLoc::recvGetLocRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad Latitude: no number, ignoring response" << std::endl;
      return -1;
    }
    if (latitude_p90 >= 180*3600) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetLoc::recvGetLocRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad Latitude: Latitude+90Deg>=" << (180*3600) << ", ignoring response" << std::endl;
      return -1;
    }

    switch (telescope.recv_buf[13]) {
      case '0':
          // no GPS module or error
        break;
      case '1':
          // GPS module ok but no GPS data (yet)
        break;
      case '2':
          // GPS data ok
        break;
      default:
        std::cout << PrintTime() << " "
                     "TelescopeIOptron::CommandGetLoc::recvGetLocRsp: \""
                  << std::string(telescope.recv_buf,telescope.recv_used)
                  << "\": bad GPS status '"
                  << telescope.recv_buf[13] << '\'' << std::endl;
        return -1;
    }
    switch (telescope.recv_buf[14]) {
      case '0':
          // stopped
        break;
      case '1':
          // tracking without PEC
        break;
      case '2':
          // slewing
        break;
      case '3':
          // auto-guiding
        break;
      case '4':
          // meridian flipping
        break;
      case '5':
          // tracking with PEC
        break;
      case '6':
          // parked
        break;
      case '7':
          // stopped at zero position
        break;
      default:
        std::cout << PrintTime() << " "
                     "TelescopeIOptron::CommandGetLoc::recvGetLocRsp: \""
                  << std::string(telescope.recv_buf,telescope.recv_used)
                  << "\": bad system status '"
                  << telescope.recv_buf[14] << '\'' << std::endl;
        return -1;
    }
    switch (telescope.recv_buf[15]) {
      case '0':
          // sideral rate
        break;
      case '1':
          // lunar rate
        break;
      case '2':
          // solar rate
        break;
      case '3':
          // King rate
        break;
      case '4':
          // custom rate
        break;
      default:
        std::cout << PrintTime() << " "
                     "TelescopeIOptron::CommandGetLoc::recvGetLocRsp: \""
                  << std::string(telescope.recv_buf,telescope.recv_used)
                  << "\": bad tracking rate '"
                  << telescope.recv_buf[15] << '\'' << std::endl;
        return -1;
    }
    telescope.locationReceived(longitude_sign ? (360*3600u-abs_longitude) : abs_longitude, latitude_p90);
    telescope.getPos();
    telescope.commandFinished();
    return 20;
  }
};

void TelescopeIOptron::getLoc(void) {
  if (!next_command_get_loc) {
    next_command_get_loc = std::make_unique<CommandGetLoc>(*this);
  }
  doSomething();
}


//------------------------------------------------------------------------



//---------------Pointing Position----------------------------------------

void TelescopeIOptron::positionReceived(const IOptronRaDec &ra_dec) {
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
                 "TelescopeIOptron::positionReceived" << ra_dec << ": "
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
                   "TelescopeIOptron::positionReceived" << ra_dec << ": "
                   "J2000(Ra:" << PrintRaInt(ra_int_j2000)
                << ",Dec:" << PrintDecInt(dec_int_j2000) << "), "
                   "H:" << PrintRaInt(hour_angle_int)
                << ", Az:" << PrintRaRad(az) << ", Alt:" << PrintDecRad(alt);
#endif

    if (alt < -2.0*(M_PI/180) || // allow observing terrestrial objects
        ForbiddenDecHA(ra_dec.getDecInt(),hour_angle_int)) {
        // GOTO hour angle 1h, keep dec
      const unsigned int goto_ra_int = W_int-longitude_int-(0x40000003u/6u);
      const int goto_dec_int = std::min(std::max(ra_dec.getDecInt(),0),(0x40000000/90)*80);

          // do not send another GOTO until ths GOTO has had a chance to finish:
      unsigned int d_ra = goto_ra_int - ra_dec.getRaInt();
      if (d_ra >= 0x80000000u) d_ra = -d_ra;
      int d_dec = goto_dec_int - ra_dec.getDecInt();
      if (d_dec < 0) d_dec = -d_dec;
      const long long int expected_duration = (std::max(d_ra,(unsigned int)d_dec)*(86400*1000000LL/1024)) >> 32;
      end_of_last_goto = now + expected_duration;

#ifdef PRINT_PERIODIC_POSITION
        std::cout << PrintRaInt(hour_angle_int)
                  << " FORBIDDEN, "
//                     "expecting " << expected_duration << "us for "
                     "GOTO("
                  << "Ra:" << PrintRaInt(goto_ra_int)
                  << ",dec:" << PrintDecInt(goto_dec_int) << ')';
#endif
      gotoPosition(IOptronRaDec(goto_ra_int,goto_dec_int));
    }
#ifdef PRINT_PERIODIC_POSITION
      std::cout << std:: endl;
#endif
  }
}


class TelescopeIOptron::CommandGetPos : public TelescopeIOptron::Command {
public:
  CommandGetPos(TelescopeIOptron &telescope) : Command(telescope) {}
  void execAsync(void) override {
    sendGetPosRqu();
  }
  void print(std::ostream &o) const override {o << "CommandGetPos";}
private:
  void sendGetPosRqu(void) {
    telescope.sendRqu(":GEC#",5,std::bind(&TelescopeIOptron::CommandGetPos::recvGetPosRsp,this));
  }
  int recvGetPosRsp(void) {
    if (telescope.recv_used < 18) return 0;
    if (telescope.recv_buf[17] != '#') {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetPos::recvGetPosRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": no '#' at end, ignoring response" << std::endl;
      return -1;
    }
    IOptronRaDec ra_dec;
    ra_dec.dec_sign = (telescope.recv_buf[0] == '-');
    if (!ra_dec.dec_sign && (telescope.recv_buf[0] != '+')) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetPos::recvGetPosRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": Dec sign must be '+' or '-', ignoring response" << std::endl;
      return -1;
    }
    if (!ParseDecimal(telescope.recv_buf+1,8,ra_dec.dec)) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetPos::recvGetPosRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad Dec: no number, ignoring response" << std::endl;
      return -1;
    }
    if (ra_dec.dec > 90*3600*100) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetPos::recvGetPosRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad Dec: abs(Dec)>" << (90*3600*100) << ", ignoring response"  << std::endl;
      return -1;
    }
    if (!ParseDecimal(telescope.recv_buf+9,8,ra_dec.ra)) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetPos::recvGetPosRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad Ra: no number, ignoring response" << std::endl;
      return -1;
    }
    if (ra_dec.ra >= 24*3600*1000) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGetPos::recvGetPosRsp: \""
                << std::string(telescope.recv_buf,telescope.recv_used)
                << "\": bad Ra: Ra>=" << (24*3600*1000) << ", ignoring response" << std::endl;
      return -1;
    }
//    std::cout << PrintTime() << " "
//                 "TelescopeIOptron::CommandGetPos::recvGetPosRsp: \""
//                << std::string(telescope.recv_buf,telescope.recv_used)
//                << "\" ok: " << ra_dec << std::endl;
    telescope.positionReceived(ra_dec);
    telescope.get_pos_deadline.expires_from_now(boost::posix_time::microseconds(100000));
    telescope.get_pos_deadline.async_wait(
      [t = &telescope](const boost::system::error_code &e) {
        if (e) return; // timer was cancelled
        t->getPos();
      });
    telescope.commandFinished();
    return 18;
  }
};

void TelescopeIOptron::getPos(void) {
  if (!next_command_get_pos) {
    next_command_get_pos = std::make_unique<CommandGetPos>(*this);
  }
  doSomething();
}

//------------------------------------------------------------------------




class TelescopeIOptron::CommandGoto : public TelescopeIOptron::Command {
public:
  CommandGoto(TelescopeIOptron &telescope) : Command(telescope) {}
  void set(const IOptronRaDec &ra_dec) {
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
    FillDecimal(buf+4,8,ra_dec.getAbsDec());
    buf[12] = '#';
    telescope.sendRqu(buf,13,std::bind(&TelescopeIOptron::CommandGoto::recvDecRsp,this));
  }
  int recvDecRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGoto" << ra_dec << "::recvDecRsp \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                   "Sd: Bad Response" << std::endl;
      return -1;
    }
    sendRaRqu();
    return 1;
  }
  void sendRaRqu(void) {
    buf[0] = ':';
    buf[1] = 'S';
    buf[2] = 'r';
    FillDecimal(buf+3,8,ra_dec.getRa());
    buf[11] = '#';
    telescope.sendRqu(buf,12,std::bind(&TelescopeIOptron::CommandGoto::recvRaRsp,this));
  }
  int recvRaRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandGoto" << ra_dec << "::recvRaRsp \""
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
    telescope.sendRqu(buf,4,std::bind(&TelescopeIOptron::CommandGoto::recvGotoRsp,this));
  }
  int recvGotoRsp(void) {
    switch (telescope.recv_buf[0]) {
      case '1':
//          std::cout << PrintTime() << " "
//                       "TelescopeIOptron::CommandGoto" << ra_dec << "::recvGotoRsp \""
//                    << std::string(telescope.recv_buf,recv_used) << "\": "
//                       "GOTO ok."
//                    << std::endl;
        break;
      case '0':
        std::cout << PrintTime() << " "
                     "TelescopeIOptron::CommandGoto" << ra_dec << "::recvGotoRsp \""
                  << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                     "GOTO not ok. The desired object is below 0 degrees altitude."
                  << std::endl;
        break;
      default:
        std::cout << PrintTime() << " "
                     "TelescopeIOptron::CommandGoto" << ra_dec << "::recvGotoRsp \""
                  << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                     "MS: Bad Response" << std::endl;
        return -1;
    }
    telescope.commandFinished();
    return 1;
  }
private:
  IOptronRaDec ra_dec;
  char buf[13];
};

void TelescopeIOptron::gotoPosition(const IOptronRaDec &ra_dec) {
  std::cout << PrintTime() << " "
               "TelescopeIOptron::gotoPosition" << ra_dec << std::endl;
  if (!next_command_goto) {
    next_command_goto = std::make_unique<CommandGoto>(*this);
  }
  next_command_goto->set(ra_dec);
  doSomething();
}

void TelescopeIOptron::gotoPosition(const unsigned int ra_int_j2000,const int dec_int_j2000) {
//    if (dec_int_j2000 < -0x40000000 || dec_int_j2000 > 0x40000000) abort();

  const Vector<double,3> v0 = PolarToRect( ra_int_j2000 * (M_PI/2147483648.0),
                                          dec_int_j2000 * (M_PI/2147483648.0) );
  const Vector<double,3> v = v0*precession_matrix;
  double ra,dec;
  RectToPolar(v,ra,dec);
  if (dec < -0.5*M_PI || dec > 0.5*M_PI) abort();
  IOptronRaDec ra_dec( (unsigned int)floor(0.5+ ra*(2147483648.0/M_PI)),
                                (int)floor(0.5+dec*(2147483648.0/M_PI)) );
  gotoPosition(ra_dec);
}

//------------------------------------------------------------------------




class TelescopeIOptron::CommandMove : public TelescopeIOptron::Command {
public:
  CommandMove(TelescopeIOptron &telescope) : Command(telescope) {}
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
//          std::cout << "CommandMove::move_deadline_l: " << e.message() << std::endl;
          if (e) return; // timer was cancelled
          std::cout << "CommandMove::move_deadline_l: stopping movement" << std::endl;
          t->move(0,0,0);
        });
    }
    if (telescope.curr_horz == horz && telescope.curr_vert == vert) {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandMove::execAsync: already stopped, nothing to do." << std::endl;
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
                   "TelescopeIOptron::CommandMove::sendSpeedRqu: stop" << std::endl;
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
    telescope.sendRqu(buf,5,std::bind(&TelescopeIOptron::CommandMove::recvSpeedRsp,this));
  }
  int recvSpeedRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandMove(" << horz << ',' << vert << ")::recvSpeedRsp \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
                   "Bad Response" << std::endl;
      return -1;
    }
    telescope.curr_speed = speed;
    sendHorzRqu();
    return 1;
  }
  void sendStopRqu(void) {
    telescope.sendRqu(":q#",3,std::bind(&TelescopeIOptron::CommandMove::recvStopRsp,this));
  }
  int recvStopRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandMove(" << horz << ',' << vert << ")::recvStopRsp \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
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
      telescope.sendRqu(buf,4,std::bind(&TelescopeIOptron::CommandMove::recvHorzStopRsp,this));
    } else {
      buf[1] = 'm';
        // horz > 0: right, east
      buf[2] = (horz < 0) ? 'w' : 'e';
      telescope.sendMsg(buf,4,std::bind(&TelescopeIOptron::CommandMove::horzRquFinished,this));
    }
  }
  void horzRquFinished(void) {
    telescope.curr_horz = horz;
    sendVertRqu();
  }
  int recvHorzStopRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandMove(" << horz << ',' << vert << ")::recvHorzStopRsp \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
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
      telescope.sendRqu(buf,4,std::bind(&TelescopeIOptron::CommandMove::recvVertStopRsp,this));
    } else {
      buf[1] = 'm';
        // vert > 0: down, south
      buf[2] = (vert < 0) ? 's' : 'n'; // :mn# actually moves south
      telescope.sendMsg(buf,4,std::bind(&TelescopeIOptron::CommandMove::vertRquFinished,this));
    }
  }
  void vertRquFinished(void) {
    telescope.curr_vert = vert;
    telescope.commandFinished();
  }
  int recvVertStopRsp(void) {
    if (telescope.recv_buf[0] != '1') {
      std::cout << PrintTime() << " "
                   "TelescopeIOptron::CommandMove(" << vert << ',' << vert << ")::recvVertStopRsp \""
                << std::string(telescope.recv_buf,telescope.recv_used) << "\": "
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

void TelescopeIOptron::move(short int horz,short int vert,
          unsigned int validity_micros) {
  if (!next_command_move) {
    next_command_move = std::make_unique<CommandMove>(*this);
  }
  next_command_move->accumulate(horz,vert,validity_micros);
  next_command_goto.reset();
  doSomething();
}



//------------------------------------------------------------------------







class TelescopeIOptron::CommandGuide  : public TelescopeIOptron::Command {
public:
  CommandGuide(TelescopeIOptron &telescope) : Command(telescope),d_ra(0),d_dec(0) {}
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
    telescope.sendMsg(buf,9,std::bind(&TelescopeIOptron::CommandGuide::raRquFinished,this));
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
    telescope.sendMsg(buf,9,std::bind(&TelescopeIOptron::CommandGuide::decRquFinished,this));
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

void TelescopeIOptron::guide(int d_ra_micros,int d_dec_micros) {
  if (!next_command_guide) {
    next_command_guide = std::make_unique<CommandGuide>(*this);
  }
  next_command_guide->accumulate(d_ra_micros,d_dec_micros);
  doSomething();
}


//------------------------------------------------------------------------




void TelescopeIOptron::init(unsigned int drain_micros) {
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

void TelescopeIOptron::commandFinished(void) {
//    std::cout << "commandFinished" << std::endl;
  command_deadline.cancel();
  curr_command.reset();
  doSomething();
}

void TelescopeIOptron::doSomething(void) {
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
//                       "TelescopeIOptron::doSomething::timeout_l: " << error.message() << std::endl;
          return; // probably cancelled
        }
        std::cout << PrintTime() << " "
                     "TelescopeIOptron::doSomething::timeout_l: " << *curr_command << " timeout" << std::endl;
        serial.cancel();
        doSomething();
      });
  }
//    std::cout << "doSomething end" << std::endl;
}



