# Makefile for integrator!

IDIR = ../headers
ODIR = objects
OUTDIR = ../output

CPP = g++
CFLAGS = -g -I $(IDIR)
LIBS = 

_OBJECTS = network.o nonaffinity.o print.o integrator.o 
OBJECTS = $(patsubst %, $(ODIR)/%, $(_OBJECTS))

_INCLUDE = network.h print.h nonaffinity.h utils.h
INCLUDE = $(patsubst %, $(IDIR)/%, $(_INCLUDE))

$(ODIR)/%.o: %.cpp $(INCLUDE)
	$(CPP) -c -o $@ $< $(CFLAGS)

integrator: $(OBJECTS)
	$(CPP) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean
	
clean: 
	-rm -f $(OUTDIR)/*.txt
