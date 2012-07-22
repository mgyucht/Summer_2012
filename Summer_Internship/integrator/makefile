# Makefile for integrator!

IDIR = src
ODIR = src/objects
SDIR = src
EXECDIR = .

HOSTNAME = $(shell hostname)

CPP = g++
CPPFLAGS = -g -O1 -I $(IDIR) -Wall
LIBS = -lboost_program_options

_OBJECTS = network.o nonaffinity.o print.o integrator.o motors.o
OBJECTS = $(patsubst %, $(ODIR)/%, $(_OBJECTS))

_INCLUDE = network.h print.h nonaffinity.h utils.h motors.h
INCLUDE = $(patsubst %, $(IDIR)/%, $(_INCLUDE))

$(ODIR)/%.o: $(SDIR)/%.cpp $(IDIR)/%.h makefile
	@mkdir -p $(ODIR)
	$(CPP) -c -o $@ $< $(CFLAGS)

integrator.out: $(OBJECTS) $(INCLUDE)
	$(CPP) $(CFLAGS) $(LIBS) -o $(EXECDIR)/$@ $^

.PHONY: clean
	
clean:
ifeq ($(HOSTNAME), della3)
	rm -rf ~/data/*
	rm ~/Summer_2012/Summer_Internship/scripts/*\.log*
else
	find output/ -type f -name "*.txt" -exec rm {} \;
	find output/ -type f -name "*.tiff" -exec rm {} \;
endif
