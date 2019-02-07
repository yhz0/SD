# Makefile definitions for CPLEX 12.8
#    10 Oct. 2018
# For Ubuntu 12.04 LTS
# Written by Harsha Gangammanavar
#------------------------------------------------------------

#------------------------------------------------------------
# System Description for My Computer
#------------------------------------------------------------
SYSTEM     = x86-64_sles10_4.1
LIBFORMAT  = static_pic

#------------------------------------------------------------
# Location of CPLEX installations:
# TODO: Replace with your installation folder
#------------------------------------------------------------
CPLEXDIR = /opt/ilog/cplex128/cplex
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXINCDIR   = $(CPLEXDIR)/include

#------------------------------------------------------------
# Compiler Selection (CC is for c)
#------------------------------------------------------------
CC  = gcc -O3

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------
COPT  = -m64 -fPIC -fno-strict-aliasing 

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
CLNFLAGS  = -lcplex -lm -pthread

# ---------------------------------------------------------------------
# To execute the CPLEX examples.
# ---------------------------------------------------------------------
CFLAGS  = $(COPT) -I$(CPLEXINCDIR) -L$(CPLEXLIBDIR)

#------------------------------------------------------------
exe = twoSD
src = $(wildcard *.c)
dep = $(wildcard *.h)
obj = $(src:.c=.o)

#SOURCES = twoSD.c algo.c cuts.c master.c prob.c rvgen.c smps.c solver.c subprob.c utility.c \
#compromise.c evaluate.c optimal.c randCost.c setup.c soln.c stocUpdate.c
# dep =prob.h smps.h solver.h stoc.h twoSD.h utils.h
#OBJECTS = twoSD.o algo.o cuts.o master.o prob.o rvgen.o smps.o solver.o subprob.o utility.o \
#compromise.o evaluate.o optimal.o randCost.o setup.o soln.o stocUpdate.o

# ------------------------------------------------------------
%.o: %.c $(dep)
	$(CC) -c -o $@ $< $(CFLAGS) $(CLNFLAGS)

$(exe): $(obj)
	$(CC) -o $(exe) $(obj) $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) $(exe)

# make definitions:
#%.o : %.c
#	$(CC) -c -o $@ $< $(CFLAGS)

#$(EXECUTABLE) : $(OBJECTS)
#	$(CC) $(OBJECTS) -o $(EXECUTABLE) $(CFLAGS) $(CLNFLAGS)

#$(OBJECTS) : $(HEADERS)
#	$(CC) -c $(SOURCES) $(CFLAGS)

#clean :
#	rm -rf sd *.o *.out *.txt *~ # *.class
#	rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp *.dat Summary *.sol
#	rm -rf Batch_dual Master_Dual Batch_Obj Master_Obj Batch_x incumb
