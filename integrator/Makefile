# Makefile for integrator!

IDIR = src
ODIR = src
SDIR = src
EXECDIR = .

HOSTNAME = $(shell hostname)
CPPFLAGS = -I $(IDIR) 

ifeq ($(HOSTNAME), della3)
	CPPFLAGS += -DDELLA3
endif

CPP = g++
LIBS = -lboost_program_options

_SOURCES = integrator.cpp network.cpp nonaffinity.cpp print.cpp motors.cpp \
	   options.cpp
SOURCES = $(patsubst %, $(SDIR)/%, $(_SOURCES))

_OBJECTS = $(_SOURCES:.cpp=.o)
OBJECTS = $(patsubst %, $(ODIR)/%, $(_OBJECTS))

_INCLUDE = network.h print.h nonaffinity.h utils.h motors.h options.h
INCLUDE = $(patsubst %, $(IDIR)/%, $(_INCLUDE))

# -------------------------------------------------------------------------#

all: integrator.out

debug: integrator-noopt.out
debug: CPPFLAGS += -DDEBUG -g

%.o: %.cpp %.h makefile
	$(CPP) -c -o $@ $< $(CPPFLAGS)

integrator.out: CPPFLAGS += -O2
integrator.out: $(OBJECTS) $(INCLUDE)
	$(CPP) $(CPPFLAGS) $(LIBS) -o $(EXECDIR)/$@ $^

integrator-noopt.out: $(OBJECTS) $(INCLUDE)
	$(CPP) $(CPPFLAGS) $(LIBS) -o $(EXECDIR)/$@ $^

# -------------------------------------------------------------------------#

.PHONY: clean clena

clena:
clean:
ifeq ($(HOSTNAME), della3)
	-rm ~/Summer_2012/Summer_Internship/scripts/*log*
	-rm ~/Summer_2012/Summer_Internship/scripts/new/*log*
else
	find output/ -type f -name "*.txt" -exec rm {} \;
	find output/ -type f -name "*.tiff" -exec rm {} \;
endif
