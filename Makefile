# Makefile definitions for CPLEX 12.8
#    10 Oct. 2018
# For Ubuntu 12.04 LTS
# Written by Harsha Gangammanavar
#------------------------------------------------------------

#------------------------------------------------------------
# TODO: Replace with appropriate folder names and system descriptions
# Location of CPLEX installations and system descriptions
CPLEXDIR = /opt/ibm/ILOG/CPLEX_Studio1210/cplex
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXINCDIR   = $(CPLEXDIR)/include
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
SUBDIRS := ./src

#------------------------------------------------------------
# Compiler, Linker Selections and Definitions (CC is for c)
OBJ_SRCS := 
ASM_SRCS := 
C_SRCS := $(wildcard $(SUBDIRS)/*.c)
O_SRCS := 
S_UPPER_SRCS := 
EXECUTABLES := 
OBJS := $(C_SRCS:.c=.o)
C_DEPS := $(C_SRCS:.c=.d)

CC  := gcc-5
LIBS := -lilocplex -lcplex -lpthread -lm -ldl
INCLUDES := -I$(CPLEXINCDIR) -I$(SUBDIRS)
CFLAGS := 

#------------------------------------------------------------
# All Target
all: twoSD

# Tool invocations
twoSD: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Cross GCC Linker'
	$(CC) $(CFLAGS) -L$(CPLEXLIBDIR) -o "twoSD" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

library: $(OBJS) $(USER_OBJS)
	$(CC) -fPIC -shared -L$(CPLEXLIBDIR) -o "libtwosd.so" $(OBJS) $(USER_OBJS) $(LIBS)

%.o: %.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	$(CC) $(CFLAGS) $(INCLUDES) -O0 -g3 -Wall -c -fmessage-length=0 -fPIC -fexceptions -m64 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" -c "$<"
	@echo 'Finished building: $<'
	@echo ' '


# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) twoSD
	-@echo ' '

.PHONY: all clean dependents
