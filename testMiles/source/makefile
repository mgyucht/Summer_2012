# Makefile for integrator!

IDIR = ../headers
ODIR = objects

CPP = g++
CFLAGS = -g -I $(IDIR)
LIBS = 

_OBJECTS = program.o network.o nonaffinity.o print.o
OBJECTS = $(patsubst %, $(ODIR)/%, $(_OBJECTS))

_INCLUDE = network.h print.h nonaffinity.h utils.h
INCLUDE = $(patsubst %, $(IDIR)/%, $(_INCLUDE))

$(ODIR)/%.o: %.cpp $(INCLUDE)
	$(CPP) -c -o $@ $< $(CFLAGS)

program: $(OBJECTS)
	$(CPP) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean
	
clean: 
	cd ../output
	-rm -f *.txt
	cd ../source
