BOOST_INCLUDE = /opt/homebrew/include
BOOST_LIB = /opt/homebrew/lib/

CC=g++
CFLAGS=-c -Wall -I./libs_gamb -I./libs_DockingLib -I$(BOOST_INCLUDE) -O2 --std=c++11 # -fexpensive-optimizations -ffast-math

SOURCES_MAIN = $(wildcard *.cc)
SOURCES_GAMB = $(wildcard libs_gamb/*.cc)
SOURCES_DOCKLIB = $(wildcard libs_DockingLib/*.cc)

OBJECTS_MAIN = $(SOURCES_MAIN:.cc=.o)
OBJECTS_GAMB = $(SOURCES_GAMB:.cc=.o)
OBJECTS_DOCKLIB = $(SOURCES_DOCKLIB:.cc=.o)

all: MainFast

MainFast: libgamb.a libdocklib.a $(OBJECTS_MAIN)
	$(CC) $(OBJECTS_MAIN) -L. -L$(BOOST_LIB) -lgamb -ldocklib -lboost_program_options -o CombinatorialAssembler.out 

%.o: %.cc
	$(CC) $(CFLAGS) $< -o $@

libgamb.a: $(OBJECTS_GAMB)
	ar rcs libgamb.a $(OBJECTS_GAMB)

libdocklib.a: $(OBJECTS_DOCKLIB) libgamb.a
	ar rcs libdocklib.a $(OBJECTS_DOCKLIB) $(OBJECTS_GAMB)

clean:
	rm -f *.o *.a CombinatorialAssembler.out libs_gamb/*.o libs_DockingLib/*.o

