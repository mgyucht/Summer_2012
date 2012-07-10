# Makefile for integrator!

IDIR = headers
ODIR = source/objects
SDIR = source
EXECDIR = .

CPP = g++
CFLAGS = -g -I $(IDIR) -Wall
LIBS = 

_OBJECTS = network.o nonaffinity.o print.o integrator.o 
OBJECTS = $(patsubst %, $(ODIR)/%, $(_OBJECTS))

_INCLUDE = network.h print.h nonaffinity.h utils.h compileinfo.h
INCLUDE = $(patsubst %, $(IDIR)/%, $(_INCLUDE))

$(ODIR)/%.o: $(SDIR)/%.cpp $(INCLUDE)
	$(CPP) -c -o $@ $< $(CFLAGS)

integrator.out: $(OBJECTS)
	$(CPP) -o $(EXECDIR)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean
	
clean: 
	rm -rf output/*
