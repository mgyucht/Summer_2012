# clean.mk
# --------
#  
#  The  clean  target  for   make   clears   output(/logs)   from   executing
#  integrator.out. It is different between my personal computer and on della,
#  so I include this file in the makefile and allow the makefile to  be  part
#  of the git repository but not this file.

clean:
	rm -rf ~/data/*
	rm ~/Summer_2012/Summer_Internship/scripts/*\.log*
