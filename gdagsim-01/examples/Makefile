# GDAGsim 0.1
# (C) 2000, Darren J Wilkinson

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
# ANSI C compiler
CC=gcc

# You may also want to tweak these
CFLAGS=-I$(GSL_INC) -I$(GDAG_INC) -ansi -Wall -pedantic -g
LDFLAGS=-lgsl -lgslblas -L$(LIBPATH) -lgdag

# **********************************************************************
# Implicit rules will do the rest.
# eg. "make example_dlm" will build executable "example_dlm" from
#     source file "example_dlm.c".
# **********************************************************************

clean:
	rm -f *~ *.o *.out core *.tab *.ind example_dlm example_twa



# eof 

