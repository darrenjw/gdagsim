# GDAGsim 0.2
# (C) 2000, Darren J Wilkinson

# Makefile for building GDAGsim library and test code


# **********************************************************************
# Edit these to suit your system
# **********************************************************************


# where you want the library installing
LIBPATH=${HOME}/lib
# where you want header installing
INCPATH=${HOME}/lib/include
# location of GSL headers
GSL_INC=/usr/include/gsl
# location of Meschach headers
MES_INC=/usr/include/meschach
# ANSI C compiler
CC=gcc

# You may also want to tweak these
CFLAGS=-I$(GSL_INC) -I$(MES_INC) -I. -Wall -pedantic -g -c
# CFLAGS=-I$(GSL_INC) -I$(MES_INC) -I. -ansi -Wall -pedantic -g -c

# GSl and BLAS libraries
LDFLAGS=-lgsl -lgslcblas -lmeschach
# linking in the new library
LDEXTRA=-L $(LIBPATH) -lgdag


# **********************************************************************
# You shouldn't need to edit anything below here
# **********************************************************************


# object files comprising the main library
OBJS=gdag_main.o gdag_sim.o gdag_utils.o gdag_sparse.o gdag_meschach.o
# C source files for the main library
SRC=gdag_main.c gdag_sim.c gdag_utils.c gdag_sparse.c gdag_meschach.c
# all other C source files
OCSRC=gdag_test.c gdag_sparse_test.c
# other source files
OSRC=Makefile gdag.h gdag_meschach.h
# all source files
SOURCE=$(OSRC) $(SRC) $(OCSRC)
# all exectutables
EXEC=gdag_test gdag_libtest gdag_sparse_test

# run test code
test: gdag_test gdag_sparse_test
	./gdag_sparse_test
	./gdag_test

libtest: gdag_libtest
	./gdag_libtest


# this is how to build the test code
gdag_test: gdag_test.o $(OBJS)
	$(CC) $(LDFLAGS) gdag_test.o $(OBJS) -o gdag_test

gdag_libtest: gdag_test.o
	$(CC) $(LDFLAGS) $(LDEXTRA) gdag_test.o -o gdag_libtest

gdag_sparse_test: gdag_sparse_test.o gdag_sparse.o
	$(CC) $(LDFLAGS) gdag_sparse_test.o gdag_sparse.o -o gdag_sparse_test



# build the shared library
libgdag.so: $(OBJS)
	$(CC) -shared -Wl,-soname,libgdag.so -o libgdag.so $(OBJS)


# general compilation of object files
%.o: %.c
	$(CC) $(CFLAGS) $< -o $@


# misc housekeeping functions
install: libgdag.so
	rm -f $(LIBPATH)/libgdag.so
	cp libgdag.so $(LIBPATH)/libgdag.so.2
	ln -s $(LIBPATH)/libgdag.so.2 $(LIBPATH)/libgdag.so
	cp gdag.h gdag_meschach.h $(INCPATH)

clean:
	rm -f *~ *.o *.so core $(EXEC)

print:
	a2ps $(SOURCE)

edit:
	gnuclient $(SOURCE) &


# eof 

