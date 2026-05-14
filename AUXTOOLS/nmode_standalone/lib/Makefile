#
#************************************************************************
#                             AMBER                                    **
#                                                                      **
#                Copyright (c) 1986, 1991, 1995                        **
#             Regents of the University of California                  **
#                       All Rights Reserved.                           ** 
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#
#
SHELL=/bin/sh
VPATH=
PREC=DPREC
INCDIR=../include
include ../config.h

OBJ=  \
sys.o wallclock.o veclib.o amopen.o mexit.o nxtsec.o rgroup.o \
rfree.o matinv.o

#  default rules for Fortran and C compilation:

.F.o:   $<
	@echo "[LIB]  FC $<"
	$(VB)$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

.c.o:
	@echo "[LIB]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

#  special rules for free-format routines:

random.o:   random.F90
	@echo "[LIB]  FC $<"
	$(VB)$(FC) $(FPPFLAGS) -c $(FREEFORMAT_FLAG) $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

veclib.o:   veclib.F90
	@echo "[LIB]  FC $<"
	$(VB)$(FC) $(FPPFLAGS) -c $(FREEFORMAT_FLAG) $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

#
all:
	@echo "(lib: make specific .o files, driven by other Makefiles)"

new2oldparm$(SFX):  new2oldparm.o nxtsec.o
	@echo "[LIB]  FC $<"
	$(VB)$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o new2oldparm$(SFX) new2oldparm.o \
		nxtsec.o $(LDFLAGS) $(AMBERLDFLAGS)


sys.a:   sys.o wallclock.o veclib.o
	@echo "[LIB]  AR $@"
	$(VB)$(AR) sys.a sys.o wallclock.o veclib.o > /dev/null
	$(VB)$(RANLIB) sys.a

sys2.a:  $(OBJ)
	@echo "[LIB]  AR $@"
	$(VB)$(AR) sys2.a $(OBJ)  > /dev/null
	$(VB)$(RANLIB) sys2.a

install: new2oldparm$(SFX)
	$(VB)$(MV) new2oldparm$(SFX) $(BINDIR)

clean:
	/bin/rm -f malloc_test new2oldparm$(SFX) *.o  *.a  *.mod

uninstall:
	/bin/rm -f $(BINDIR)/new2oldparm$(SFX)

