CC = g++
CFLAGS = -c -Wall -std=c++11
CPPFLAGS = -DVLEER -DBC_P1 -DBC_P2 #-DMHD -DSHEAR -DSELF_GRAVITY
SOURCES = zeus_main.cpp \
					constants.cpp \
					grids.cpp \
					source.cpp \
					timestep.cpp \
					bcs.cpp \
					transport.cpp
OBJECTS = ${SOURCES:.cpp=.o}
#EXEDIR = ../bin/
EXEC = main

all: main

main: ${OBJECTS}
	${CC} ${OBJECTS} -o ${EXEC}

zeus_main.o: zeus_main.cpp
	${CC} ${CFLAGS} ${CPPFLAGS} $<

constants.o: constants.cpp
	${CC} ${CFLAGS} ${CPPFLAGS} $<

grids.o: grids.cpp
	${CC} ${CFLAGS} ${CPPFLAGS} $<

source.o: source.cpp
	${CC} ${CFLAGS} ${CPPFLAGS} $<

timestep.o: timestep.cpp
	${CC} ${CFLAGS} ${CPPFLAGS} $<

bcs.o: bcs.cpp
	${CC} ${CFLAGS} ${CPPFLAGS} $<

transport.o: transport.cpp
	${CC} ${CFLAGS} ${CPPFLAGS} $<

.PHONY: clean
clean:
	rm *.o

#all: grid_specs.cpp subroutines.cpp zeus_main.cpp
#	g++ -c grid_specs.cpp
#	g++ -c subroutines.cpp
#	g++ -c zeus_main.cpp
#	g++ *.o -o zeus_main

#subs: subroutines.cpp
#	g++ -c subroutines.cpp

#grid: grids.cpp
#	g++ -c grids.cpp

#main: zeus_main.cpp init.cpp grids.cpp
#	g++ zeus_main.cpp grids.cpp init.h grids.h -std=c++11 -o zeus_main

