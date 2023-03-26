BOOST_INCLUDE = /opt/homebrew/include
BOOST_LIB = /opt/homebrew/lib/

CC=g++
# use -Wno-deprecated-declarations to suppress warnings from boost
# CFLAGS=-c -Wall -I./libs_gamb -I./libs_DockingLib -I$(BOOST_INCLUDE) -O2 --std=c++11 # -fexpensive-optimizations -ffast-math
CFLAGS=-c -Wall -Wno-deprecated-declarations -I./libs_gamb -I./libs_DockingLib -I$(BOOST_INCLUDE) -g -O2 --std=c++17 # -fexpensive-optimizations -ffast-math
# CFLAGS=-c -Wall -I./libs_gamb -I./libs_DockingLib -I$(BOOST_INCLUDE) -O0 -g --std=c++11 # -fexpensive-optimizations -ffast-math

SOURCES_MAIN = $(wildcard *.cc)
SOURCES_GAMB = $(wildcard libs_gamb/*.cc)
SOURCES_DOCKLIB = $(wildcard libs_DockingLib/*.cc)
SOURCES_AF2TRANS = $(wildcard AF2trans/*.cc)

OBJECTS_MAIN = $(SOURCES_MAIN:.cc=.o)
OBJECTS_GAMB = $(SOURCES_GAMB:.cc=.o)
OBJECTS_DOCKLIB = $(SOURCES_DOCKLIB:.cc=.o)
OBJECTS_AF2TRANS = $(SOURCES_AF2TRANS:.cc=.o)

all: MainCombAssemble MainAf2trans

MainCombAssemble: libgamb.a libdocklib.a $(OBJECTS_MAIN)
	$(CC) $(OBJECTS_MAIN) -L. -L$(BOOST_LIB) -lgamb -ldocklib -lboost_program_options -lpthread -o CombinatorialAssembler.out 

MainAf2trans: libgamb.a libdocklib.a $(OBJECTS_AF2TRANS)
	$(CC) $(OBJECTS_AF2TRANS) -L. -L$(BOOST_LIB) -lgamb -ldocklib -lboost_program_options -lpthread -o AF2trans.out 

%.o: %.cc
	$(CC) $(CFLAGS) $< -o $@

libgamb.a: $(OBJECTS_GAMB)
	ar rcs libgamb.a $(OBJECTS_GAMB)

libdocklib.a: $(OBJECTS_DOCKLIB) libgamb.a
	ar rcs libdocklib.a $(OBJECTS_DOCKLIB) $(OBJECTS_GAMB)

clean_all:
	rm -f *.o *.a AF2trans.out CombinatorialAssembler.out AF2trans/*.o libs_gamb/*.o libs_DockingLib/*.o

clean:
	rm -f *.o AF2trans/*.o AF2trans.out CombinatorialAssembler.out

