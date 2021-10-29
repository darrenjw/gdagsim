# GDAGsim 0.3
# (C) 2000-2002, Darren J Wilkinson

# Makefile for building code using the GDAGsim library

# **********************************************************************
# Edit these to suit your system
# **********************************************************************

# where the library is installed
LIBPATH=${HOME}/lib
# location of gdag.h header file
GDAG_INC=${HOME}/lib/include
# location of GSL headers
GSL_INC=/usr/include/gsl
MES_INC=/usr/include/meschach
# ANSI C compiler
CC=gcc

# You may also want to tweak these
CFLAGS=-I$(GSL_INC) -I$(MES_INC) -I$(GDAG_INC) -Wall -pedantic -g
LDFLAGS=-lgsl -lgslcblas -lmeschach -L$(LIBPATH) -lgdag


spat_test: spat_test.o spat.o
	$(CC) $(CFLAGS) $(LDFLAGS) spat_test.o spat.o -o spat_test

sptkriging: sptkriging.o spat.o
	$(CC) $(CFLAGS) $(LDFLAGS) sptkriging.o spat.o -o sptkriging

mcmc: mcmc.o spat.o
	$(CC) $(CFLAGS) $(LDFLAGS) mcmc.o spat.o -o mcmc




clean:
	rm -f *~ *.o *.out core *.tab *.ind spat_test mcmc sptkriging



# eof 

