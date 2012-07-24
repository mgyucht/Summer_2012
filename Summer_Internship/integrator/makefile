# Makefile for integrator!

IDIR = src
ODIR = src
SDIR = src
EXECDIR = .

HOSTNAME = $(shell hostname)

CPP = g++
CPPFLAGS = -g -O2 -I $(IDIR) -Wall -Weffc++ -pedantic  \
    -pedantic-errors -Wextra  -Wall -Waggregate-return -Wcast-align \
    -Wcast-qual  -Wchar-subscripts  -Wcomment -Wconversion \
    -Wdisabled-optimization \
    -Wfloat-equal  -Wformat  -Wformat=2 \
    -Wformat-nonliteral -Wformat-security  \
    -Wformat-y2k \
    -Wimport  -Winit-self  -Winline \
    -Winvalid-pch   \
    -Wunsafe-loop-optimizations  -Wlong-long -Wmissing-braces \
    -Wmissing-field-initializers -Wmissing-format-attribute   \
    -Wmissing-include-dirs -Wmissing-noreturn \
    -Wpacked  -Wpadded -Wparentheses  -Wpointer-arith \
    -Wredundant-decls -Wreturn-type \
    -Wsequence-point  -Wshadow -Wsign-compare  -Wstack-protector \
    -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch  -Wswitch-default \
    -Wswitch-enum -Wtrigraphs  -Wuninitialized \
    -Wunknown-pragmas  -Wunreachable-code -Wunused \
    -Wunused-function  -Wunused-label  -Wunused-parameter \
    -Wunused-value  -Wunused-variable  -Wvariadic-macros \
    -Wvolatile-register-var  -Wwrite-strings # -Werror
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

.PHONY: clean
	
clean:
ifeq ($(HOSTNAME), della3)
	rm -rf ~/data/*
	rm -f ~/Summer_2012/Summer_Internship/scripts/*log*
else
	find output/ -type f -name "*.txt" -exec rm {} \;
	find output/ -type f -name "*.tiff" -exec rm {} \;
endif
