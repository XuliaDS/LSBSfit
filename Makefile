# use gcc compiler
CC   = gcc

############################
###      DIRECTORIES     ###
############################

# SDIR: directory containing the source files (.c)
SDIR   = src

# HDIR: directory containing the header files (.h)
HDIR = head

# list of flags to pass to the compilation command: Call make debug to compile in debug mode
CFLAGS_RUN      = -O3 -g -std=c99 $(I_DEF) $(I_DEBUG) $(I_IEEE) -fopenmp $(I_HEAD) $(I_EGADS) -fPIC
CFLAGS_DEBUG    = -Og -std=c99 $(I_DEF) $(I_DEBUG2) $(I_DEBUG) $(I_IEEE) -fopenmp $(I_HEAD) $(I_EGADS) -fPIC

I_DEBUG = -Wno-unused-variable -Wno-maybe-uninitialized -Wno-unused-but-set-variable 
I_DEBUG2 = -g -ggdb -gdwarf-4 -Werror=implicit-function-declaration
I_IEEE  = -frounding-math -fsignaling-nans
I_HEAD  = -I$(HDIR)
I_EGADS = -I$(ESP_ROOT)/include
I_DEF   = -DSTANDALONE

# must include the math library
LIBS_RUN   = $(L_EGADS) -ldl -lpthread -lm

L_EGADS = -L$(ESP_ROOT)/lib/ -legads
############################
###        FILES         ###
############################

### Source files ###
# source_temp will contain all files from the SDIR directory 
source_temp   =  $(wildcard $(SDIR)/*.c) 

# specify explicitely which source files should not be compiled
#src/approxML.c
excluded_sources = backLSfit.c







# SRC will contain only those files which remain
SRC    = $(filter-out $(excluded_sources),$(source_temp))

### Object files ###
# OBJ will contain the same file titles as SRC with the object extansion
OBJ    = $(SRC:.c=.o)

### Header files ###
# declare header files
excluded_headers = 

# add path to headers
header_temp = $(wildcard $(HDIR)/*.h)
HEAD        = $(filter-out $(excluded_headers),$(header_temp))

### Executable files ###
EXE         = LSBsplineFit

############################
###   RULES DEFINITION   ###
############################

# Target to make in running mode
# Call "make run" or "make" to use
.PHONY: run
run: CFLAGS = $(CFLAGS_RUN)
run: LIBS   = $(LIBS_RUN)
run: $(EXE)

# Target to make with debugging option
# Call "make debug" or "make" to use
.PHONY: debug
debug: CFLAGS = $(CFLAGS_DEBUG)
debug: LIBS   = $(LIBS_RUN)
debug: $(EXE)

# rule to make the object files
# source files and header files are prerequisites
%.o: %.c  $(HEAD)
	$(CC) $(CFLAGS) -c -o $@ $<

# rule to make the standalone executable
# object files are prerequisites
$(EXE): $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $(OBJ)  $(PYOBJ) $(LIBS)

# prevent make from doing something with some file names
.PHONY: clean clean_all

# define cleaning rule: remove EXEcutable, object files and python libraries
clean:
	rm -f $(OBJ)

clean_all:
	make clean
	rm -f $(EXE)



