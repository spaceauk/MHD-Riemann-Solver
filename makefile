SHELL = /bin/sh

# Select compiler
CMPR = gfortran 

CPPFLAGS=$(ARG)
#======================================================================
CPPFLAGS += -P -DOPENMP
#======================================================================
CPP = 

# Compiler-related flags
CFLAGS = -fdefault-real-8 
CFLAGS += -fcheck=all 
# Linking related flags
# LFLAGS = 

ifneq (,$(findstring OPENMP,$(CPPFLAGS)))
  CFLAGS += -fopenmp
  LFLAGS = -fopenmp
endif

OBJ = module.o input.o IC2Dtype.o MUSCL2D.o slopelimiter.o riemannS.o CT2D_variant.o CT2D.o main.o check_2Ddata.o mhddef.o write_timings.o 
EXEC = main.x

# Add other subroutines
OBJ += WagaanPS.o
OBJ += Cada3rd.o
OBJ += posiCada_mod.o

# Compile the subroutines
%.o:%.F90
	$(CMPR) -c $(CFLAGS) $<

# Link into an executable
main: $(OBJ)
	$(CMPR) $(OBJ) $(LFLAGS) -o $(EXEC)

clean:
	rm -f ./*.o ./*.x ./Error_File ./Output_File
