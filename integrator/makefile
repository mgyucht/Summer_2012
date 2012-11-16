# Makefile for integrator!

IDIR = src
ODIR = src
SDIR = src
EXECDIR = .

HOSTNAME = $(shell hostname)

CPP = g++
CPPFLAGS = -g -I $(IDIR) # -O2 -Wall -Weffc++ -pedantic  \
    -pedantic-errors -Wextra
LIBS = -lboost_program_options

_SOURCES = integrator.cpp network.cpp nonaffinity.cpp print.cpp motors.cpp
SOURCES = $(patsubst %, $(SDIR)/%, $(_SOURCES))

_OBJECTS = $(_SOURCES:.cpp=.o)
OBJECTS = $(patsubst %, $(ODIR)/%, $(_OBJECTS))

_INCLUDE = network.h print.h nonaffinity.h utils.h motors.h
INCLUDE = $(patsubst %, $(IDIR)/%, $(_INCLUDE))

%.o: %.cpp %.h makefile
	$(CPP) -c -o $@ $< $(CPPFLAGS)

integrator.out: $(OBJECTS) $(INCLUDE)
	$(CPP) $(CPPFLAGS) $(LIBS) -o $(EXECDIR)/$@ $^

.PHONY: clean clena

clena:
clean:
ifeq ($(HOSTNAME), della3)
	find ~/data/ -name "*.txt" -exec rm {} \;
	rm ~/Summer_2012/Summer_Internship/scripts/*log*
else
	find output/ -type f -name "*.txt" -exec rm {} \;
	find output/ -type f -name "*.tiff" -exec rm {} \;
endif
