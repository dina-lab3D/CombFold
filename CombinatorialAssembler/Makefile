CC=g++
CFLAGS=-c -Wall -I./libs_gamb -I./libs_DockingLib -O2 # -fexpensive-optimizations -ffast-math

SOURCES_MAIN = $(wildcard *.cc)
SOURCES_GAMB = $(wildcard libs_gamb/*.cc)
SOURCES_DOCKLIB = $(wildcard libs_DockingLib/*.cc)

OBJECTS_MAIN = $(SOURCES_MAIN:.cc=.o)
OBJECTS_GAMB = $(SOURCES_GAMB:.cc=.o)
OBJECTS_DOCKLIB = $(SOURCES_DOCKLIB:.cc=.o)

all: MainFast

MainFast: $(OBJECTS_MAIN) libgamb.a libdocklib.a
	$(CC) $(OBJECTS_MAIN) -L. -lgamb -ldocklib -lstdc++ -lboost_program_options -lX11 -lm -lpthread -o CombDock.Linux

%.o: %.cc
	$(CC) $(CFLAGS) $< -o $@

libgamb.a: $(OBJECTS_GAMB)
	ar rcs libgamb.a $(OBJECTS_GAMB)

libdocklib.a: $(OBJECTS_DOCKLIB) libgamb.a
	ar rcs libdocklib.a $(OBJECTS_DOCKLIB) $(OBJECTS_GAMB)

clean:
	rm -f *.o *.a CombDock.Linux libs_gamb/*.o libs_DockingLib/*.o

