# Makefile for integrator!

CC = g++
CFLAGS = -c -g
LFLAGS = 
OBJECTS = integrator.o network.o nonaffinity.o print.o

all: integrator

integrator: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o integrator
	
integrator.o: integrator.cpp network.h utils.h 
	$(CC) $(CFLAGS) integrator.cpp
	
network.o: network.cpp network.h
	$(CC) $(CFLAGS) network.cpp
	
nonaffinity.o: nonaffinity.cpp nonaffinity.h
	$(CC) $(CFLAGS) nonaffinity.cpp
	
print.o: print.cpp print.h
	$(CC) $(CFLAGS) print.cpp
	
clean: 
	-rm  -f position_data.txt energy_data.txt nonaff_data.txt
