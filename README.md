Summer\_2012.git
================

Primary Author: Miles Yucht
---------------------------

This git repository contains my work during the summer of 2012. This
includes code written in C++ to create a model consisting of a triangular
network of springs for networks of actin filaments in the cytoskeleton
of cells. Additionally, there are several papers on the topic of modeling
cell cytoskeletons and experiments measuring the mechanical properties of
cytoskeletons.

There are two simulations set up: first, the static simulation, where the
strain is unchanging during the simulation, is in the subdirectory
Summer\_Internship/testMiles. program uses the conjugate gradient algorithm to
minimize the energy of the spring network, which it can then use to calculate
the shear modulus of the network. The main file is located in program.cpp and
is compiled by the following command:

    $ g++ program.cpp -o program
    
So, to run the program, run 

    $ ./program

For help, run

    $ ./program help

The help file will tell you the correct usage to set the parameters for the
simulation. 

The second simulation is a dynamic simulation, where the strain applied to the
system is sinusoidal. This is located in the Summer\_Internship/integrator
subdirectory. For each node in the network, the net force is calculated which
corresponds to the node's motion. To compile:

    $ make integrator.out
    
To run:

    $ ./integrator.out

See the options with 

    $ ./integrator.out --help

Looking through my code, you might believe that it was written by an
unexperienced, amateur coder, and you'd be exactly right. This project has
also been a challenge in terms of learning C++ and object-oriented
programming. In each subsequent program, I strive to write code that more
closely subscribes to the paradigms of makefiles and C++.
