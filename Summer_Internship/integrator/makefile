# Makefile for integrator!

IDIR = headers
ODIR = source/objects
SDIR = source
EXECDIR = .

HOSTNAME = $(shell hostname)

LIBS = -lboost_program_options
CPP = g++
CFLAGS = -g -O2 -I $(IDIR) -Wall

_OBJECTS = network.o nonaffinity.o print.o integrator.o 
OBJECTS = $(patsubst %, $(ODIR)/%, $(_OBJECTS))

_INCLUDE = network.h print.h nonaffinity.h utils.h
INCLUDE = $(patsubst %, $(IDIR)/%, $(_INCLUDE))

$(ODIR)/%.o: $(SDIR)/%.cpp $(INCLUDE) makefile
	@mkdir -p $(ODIR)
	$(CPP) -c -o $@ $< $(CFLAGS)

integrator.out: $(OBJECTS)
	$(CPP) -o $(EXECDIR)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean
	
clean:
ifeq ($(HOSTNAME), della3)
	rm -rf ~/data/*
	rm ~/Summer_2012/Summer_Internship/scripts/*log*
else
	find output/ -type f -name "*.txt" -exec rm {} \;
	find output/ -type f -name "*.tiff" -exec rm {} \;
endif
