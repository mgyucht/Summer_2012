# Makefile for integrator!

IDIR = headers
ODIR = source/objects
SDIR = source
OUTDIR = output
EXECDIR = .

CPP = g++
CFLAGS = -g -I $(IDIR)
LIBS = 

_OBJECTS = network.o nonaffinity.o print.o integrator.o 
OBJECTS = $(patsubst %, $(ODIR)/%, $(_OBJECTS))

_INCLUDE = network.h print.h nonaffinity.h utils.h
INCLUDE = $(patsubst %, $(IDIR)/%, $(_INCLUDE))

$(ODIR)/%.o: $(SDIR)/%.cpp $(INCLUDE)
	$(CPP) -c -o $@ $< $(CFLAGS)

integrator: $(OBJECTS)
	$(CPP) -o $(EXECDIR)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean
	
clean: 
	-rm -f $(OUTDIR)/*.txt
