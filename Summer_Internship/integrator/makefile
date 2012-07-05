# Makefile for integrator!

IDIR = headers
ODIR = source/objects
SDIR = source
EXECDIR = .
_OUTDIR = output /scratch/gpfs/myucht
OUTDIR = $(patsubst %, %/*.txt, $(_OUTDIR))

CPP = g++
CFLAGS = -g -I $(IDIR) -Wall
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
	rm -rf ~/data/*
	find ../scripts -type f -name "*.log" -delete
	find ../scripts -type f -name "*.log-*" -delete
