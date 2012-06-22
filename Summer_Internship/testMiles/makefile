# makefile
# --------
#  Makefile for spring network project.

CC = g++
PROGRAM_H_FILES = program.cpp bracket.h dbrent.h debug.h dF1dim.h dlinmin.h frprmn.h funcd.h utils.h
DATAFILES = position_data.txt energy_data.txt nonaff_data.txt

all: program

program: program.o nonaffinity.o
	$(CC) program.o nonaffinity.o -o program
	
batch: batch.o nonaffinity.o
	$(CC) batch.o nonaffinity.o -o batch
	
program.o: program.cpp $(PROGRAM_H_FILES)
	$(CC) -c program.cpp
	
nonaffinity.o: nonaffinity.cpp nonaffinity.h
	$(CC) -c nonaffinity.cpp
	
batch.o: batch.cpp $(PROGRAM_H_FILES)
	$(CC) -c batch.cpp
	
clean:
	-rm $(DATAFILES)
