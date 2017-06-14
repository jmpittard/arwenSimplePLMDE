# Makefile for ARWEN (Astrophysical Research With Enhanced Numerics) code
#
# Compile: make
#
# To make use of the gdb (lldb) debugger type: gdb (lldb) arwen, run
#
# Profiling:
# The GNU profiler (pg) doesn't work with OSX.
# The google profiler does work, but requires a bit of setting up:
#  1) Download google-perftools (sudo port install google-perftools)
#  2) Include <gperftools/profiler.h>
#  3) Add ProfilerStart("<filename>.prof"); and ProfilerStop(); in the code.
#  4) Add -lprofiler to the LIB statement, and -Wl,-no_pie to LDFLAGS (this is
#     because Lion introduced Address space layout randomization (a good thing!)
#     so pprof cannot decode the random addresses and locate the appropriate function - it then
#     shows function addresses rather than names. See https://github.com/gperftools/gperftools/issues/363
#  5) Make and run the code as normal. This creates the <filename>.prof file.
#  6) To see a text output: pprof --text ./arwen <filename>.prof 
#     Currently I do not have dot installed, so the graphical output doesn't work.
#     See http://goog-perftools.sourceforge.net/doc/cpu_profiler.html for more info.
#
# Julian Pittard: 17.11.11
# Last modified:  17.11.11
#
# Libraries and CFLAGS (O3 is optimisation, pg is the GNU profiler, lprofiler is the google profiler).....
#-----------------------------------------------------------------------------
# MACBOOK
LIB = -lcfitsio
CFLAGS = -O3 
#-----------------------------------------------------------------------------
# LINUX WORKSTATION
#LIB = -lcfitsio
#CFLAGS = -O3 -I/usr/include/cfitsio
#-----------------------------------------------------------------------------

#Compilers
CCOMP = g++

#Linking options (g = debugging, pg = GNU profiler)
LDFLAGS = -O3 -fast
#LDFLAGS = -g
#LDFLAGS = -O3
#LDFLAGS = -pg -g

#OBJECTS
OBJS=arwen.o computeloop.o dtcon.o setup.o update.o writeData.o

#Targets

compile: $(OBJS)                           #This shows the dependencies
	$(CCOMP) $(LDFLAGS) -o arwen $(OBJS) $(LIB)  #This is the linker 

#Compile (use suffix rule to describe how to make a .o file from a .c file)
.cpp.o:
	$(CCOMP) $(CFLAGS) -c $<

#Dependencies (include any necessary .h files at end of lines)


#clean
clean:
	rm -f arwen $(OBJS)
