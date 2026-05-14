#  Amber configuration file.

###############################################################################

# (1)  Location of the installation

BASEDIR=./
BINDIR=./
LIBDIR=./

SHELL=/bin/sh
INSTALLTYPE=serial

CC=gcc
CFLAGS=-fPIC  -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE   
CNOOPTFLAGS=
COPTFLAGS= -O3 -mtune=native 
WARNFLAGS=-Wall -Wno-unused-function

FP_FLAGS=
LDFLAGS= 

LEX=   flex
YACC=   bison -y
AR=    ar rv
M4=    m4
RANLIB=ranlib
VB=@

#  Set the C-preprocessor.  Code for a small preprocessor is in
#    ucpp-1.3;  it gets installed as $(BINDIR)/ucpp;
CPP=ucpp -l

#  Information about Fortran compilation:
FC=gfortran
FFLAGS= -fdefault-integer-8 -fPIC  -I$(INCDIR) -lblas -llapack
FNOOPTFLAGS=-O0
FOPTFLAGS=-O3 -mtune=native 
FREEFORMAT_FLAG=-ffree-form
LM=-lm
FPP=cpp -traditional -P
FPPFLAGS= 
FCREAL8=-fdefault-real-8
NOFORTRANMAIN=-lgfortran -w
FWARNFLAGS=-Wall -Wno-unused-function


