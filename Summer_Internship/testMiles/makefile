# makefile
# --------
#  Makefile for spring network project.

CC = g++
CFLAGS = -O2
PROGRAM_H_FILES = bracket.h dbrent.h debug.h dF1dim.h dlinmin.h frprmn.h funcd.h utils.h
INTEGRATE_H_FILES = debug.h utils.h funcd.h
DATAFILES = position_data.txt energy_data.txt nonaff_data.txt

all: program

program: program.o nonaffinity.o
	$(CC) program.o nonaffinity.o -o program
	
batch: batch.o nonaffinity.o
	$(CC) batch.o nonaffinity.o -o batch
	
integrate: integrate.o nonaffinity.o
	$(CC) integrate.o nonaffinity.o -o integrate
	
program.o: program.cpp $(PROGRAM_H_FILES)
	$(CC) $(CFLAGS) -c program.cpp
	
nonaffinity.o: nonaffinity.cpp nonaffinity.h
	$(CC) $(CFLAGS) -c nonaffinity.cpp
	
batch.o: batch.cpp $(PROGRAM_H_FILES)
	$(CC) $(CFLAGS) -c batch.cpp
	
integrate.o: integrate.cpp $(INTEGRATE_H_FILES)
	$(CC) $(CFLAGS) -c integrate.cpp
	
clean:
	-rm -f $(DATAFILES)
