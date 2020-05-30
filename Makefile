TelescopeServer2_SOURCES = TelescopeServer2.C PrintRaDec.C Time.C Telescope.C TelescopeDummy.C TelescopeIOptron.C TelescopeLX200.C

TelescopeServer2_OBJECTS = $(patsubst %.c,%.o,$(patsubst %.C,%.o,$(TelescopeServer2_SOURCES)))

TelescopeClient2_SOURCES = TelescopeClient2.C

TelescopeClient2_OBJECTS = $(patsubst %.c,%.o,$(patsubst %.C,%.o,$(TelescopeClient2_SOURCES)))

TelescopeJoystickClient_SOURCES = TelescopeJoystickClient.C PrintRaDec.C Time.C 

TelescopeJoystickClient_OBJECTS = $(patsubst %.c,%.o,$(patsubst %.C,%.o,$(TelescopeJoystickClient_SOURCES)))

COMPILEFLAGS = -std=c++14 -g -O0 -Wall -Wextra -Wno-unused-parameter -Wno-unused-function

.SUFFIXES: .C .c .o .obj

.c.o:
	gcc $(COMPILEFLAGS) $< -MMD -MF $*.d -c -o $@

.C.o:
	g++ $(COMPILEFLAGS) $< -MMD -MF $*.d -c -o $@

all: TelescopeServer2 TelescopeClient2 TelescopeJoystickClient

echo:
	echo $(TelescopeServer2_OBJECTS)

TelescopeServer2: $(TelescopeServer2_OBJECTS)
	g++ $(TelescopeServer2_OBJECTS) -o $@ -lpthread

TelescopeClient2: $(TelescopeClient2_OBJECTS)
	g++ $(TelescopeClient2_OBJECTS) -o $@ -lpthread

TelescopeJoystickClient: $(TelescopeJoystickClient_OBJECTS)
	g++ $(TelescopeJoystickClient_OBJECTS) -o $@ -lpthread

clean:
	rm -f \
  $(TelescopeServer2_OBJECTS) $(TelescopeServer2_OBJECTS:.o=.d) \
  $(TelescopeClient2_OBJECTS) $(TelescopeClient2_OBJECTS:.o=.d) \
  $(TelescopeJoystickClient_OBJECTS) $(TelescopeJoystickClient_OBJECTS:.o=.d)

-include $(TelescopeServer2_OBJECTS:.o=.d)
-include $(TelescopeClient2_OBJECTS:.o=.d)
-include $(TelescopeJoystickClient_OBJECTS:.o=.d)
