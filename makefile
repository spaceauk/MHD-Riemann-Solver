SHELL = /bin/sh

# Select compiler
CMPR = gfortran 

CPPFLAGS=$(ARG)
#======================================================================
CPPFLAGS += -P -DOPENMP
#======================================================================
CPP = /lib/cpp -traditional

SRCDIR = .srcmake
#OBJBIN = bin

# Compiler-related flags
CFLAGS = -fdefault-real-8 
CFLAGS += -fcheck=all # & idev to help with debugging <-https://stackoverflow.com/questions/23130045/fortran-90-segmentation-fault-invalid-memory-reference-with-scalable-3d-ar
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


#OBJSRC=$(addprefix $(OBJBIN)/, $(OBJEUL) $(OBJFFT) )
#$(OBJSRC): module.F90

ifneq (exists, $(shell if [ -d $(SRCDIR) ]; then echo exists; fi))
        @echo "Error: $(SRCDIR) folders have not been created"
	@echo "Use the command 'make srcmake'"; exit 1
endif	

# Compile the subroutines
%.o:%.F90
	$(CPP) $(CPPFLAGS) $< $(SRCDIR)/$<
	$(CMPR) -c $(CFLAGS) $(SRCDIR)/$<

# Link into an executable
main: $(OBJ)
	$(CMPR) $(OBJ) $(LFLAGS) -o $(EXEC)

srcmake:
	@echo deleting $(SRCDIR)
	rm -rf $(SRCDIR)
	@echo constructing $(SRCDIR)
	
	mkdir -v -p $(SRCDIR)
#	mkdir -v -p $(OBJBIN)
	@echo "srcmake folders generated: code can now be compiled"

clean:
	rm -f ./*.o ./*.x ./Error_File ./Output_File

cleanest:
	- /bin/rm *.o *.mod *~ \
	@echo files cleaned
	rm -rf $(SRCDIR)
	@echo removed $(SRCDIR) 
