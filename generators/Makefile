CC = g++
CFLAGS = -c -Wall -std=c++11
# Make sure all flags are the same as the main makefile
CPPFLAGS = -DVLEER -DBC_P1 -DBC_P2 #-DMHD -DSHEAR -DSELF_GRAVITY

#make target
%: %.cpp
	${CC} ../constants.o ../grids.o ../bcs.o ../io.o ../timestep.o $@.cpp -o $@

