# Makefile for integrator!

CC = g++
CFLAGS = -c -g
LFLAGS = 

all: integrator

integrator: integrator.o network.o nonaffinity.o
	$(CC) $(LFLAGS) integrator.o network.o nonaffinity.o -o integrator
	
integrator.o: integrator.cpp debug.h network.h utils.h
	$(CC) $(CFLAGS) integrator.cpp
	
network.o: network.cpp network.h
	$(CC) $(CFLAGS) network.cpp
	
nonaffinity.o: nonaffinity.cpp nonaffinity.h
	$(CC) $(CFLAGS) nonaffinity.cpp
	
clean: 
	-rm  -f position_data.txt energy_data.txt nonaff_data.txt
