# Makefile for integrator!

IDIR = headers
ODIR = source/objects
SDIR = source
EXECDIR = .

HOSTNAME = $(shell hostname)

LIBS = 
CPP = g++
CFLAGS = -O2 -I $(IDIR) -Wall
ifeq ($(HOSTNAME), miles)
CFLAGS += -D DELLA=0 -D HOME_COMPUTER=1
else 
CFLAGS += -D DELLA=1 -D HOME_COMPUTER=0
endif


_OBJECTS = network.o nonaffinity.o print.o integrator.o 
OBJECTS = $(patsubst %, $(ODIR)/%, $(_OBJECTS))

_INCLUDE = network.h print.h nonaffinity.h utils.h compileinfo.h
INCLUDE = $(patsubst %, $(IDIR)/%, $(_INCLUDE))

$(ODIR)/%.o: $(SDIR)/%.cpp $(INCLUDE)
	@mkdir -p $(ODIR)
	$(CPP) -c -o $@ $< $(CFLAGS)

integrator.out: $(OBJECTS)
	$(CPP) -o $(EXECDIR)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean
	
clean:
ifeq ($(HOSTNAME), miles)
	rm -rf output/*
else
	rm -rf ~/data/*
	rm ~/Summer_2012/Summer_Internship/scripts/*\.log*
endif
