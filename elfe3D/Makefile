# choose compiler
FCC = ifx
FCC = flang
FCC = nvfortran
FCC = gfortran

# determine operating system
OS = $(shell uname)
# determine architecture
ARCH = $(shell uname -m)
# determine hostname
HOST = $(shell hostname)

# object files created during compilation
OBJ =   mod_constant.o \
	mod_types_basic.o \
	mod_util.o \
	mod_vector_products.o \
	mod_define_model.o \
	mod_read_mesh.o \
	mod_calculate_matrices.o \
	mod_model_parameters.o \
	mod_interp_functions.o \
	mod_calculate_local_left.o \
	mod_calculate_global_source.o \
	mod_BC.o \
	mod_sparse_matrix_operations.o \
	mod_solvers.o \
	mod_calculate_tf.o \
	mod_error_estimates.o \
	mod_tetgen_operations.o \
	elfe3d.o


# link MUMPS libaries here!
OBJUTIL = -L$(MUMPS_ROOT)/lib -lpthread -lm -lzmumps -lmumps_common -lesmumps -lpord -lmpiseq -L$(SCOTCH_ROOT)/lib -L$(METIS_ROOT)/lib -lscotch -lscotcherr -lmetis -lopenblas

##-------------------------------------------------------------------------------------
# Section with compiler flags for shared library and compiler options
##
## for gfortran
##
ifeq ($(findstring gfortran,${FCC}),gfortran) 

# compiler flags
FFLAGS = -Ofast -march=native -mtune=native -mfma -malign-data=cacheline -ffree-form -ftree-vectorize -ftree-vectorizer-verbose=2 -fopt-info-vec -fopenmp 
FFLAGS += -pg 
#
# for warnings you can include here: 
# -pedantic -Wall -Wopenmp-simd -Wvector-operation-performance

# compile-time directives (for preprocessor etc.)
#COMPFLAGS = -x f95-cpp-input

#Intel OpenAPI compiler
else ifeq ($(findstring ifx,${FCC}),ifx)

FFLAGS = -fast -march=x86-64-v3 -free -fopenmp -shared -fPIC
FFLAGS = -O3 -free -qopenmp -g -init=zero -mcmodel=large
#FFLAGS = -O0 -free -qopenmp -g -init=zero 

#Intel older ifort compiler
else ifeq ($(findstring ifort,${FCC}),ifort)

FFLAGS = -fast -march=core-avx2 -free -qopenmp -shared -fPIC
FFLAGS = -O3 -march=core-avx2 -free -qopenmp -fPIC

#AMD aocc compiler
else ifeq ($(findstring flang,${FCC}),flang)

FFLAGS = -Ofast -ffree-form -fopenmp 

#Nvidia nvhpc compiler
else ifeq ($(findstring nvfortran,${FCC}),flang)

FFLAGS = -Ofast -Mfreeform -mp 

else 
endif

FFLAGS += -I$(MUMPS_ROOT)/include

##-------------------------------------------------------------------------------------
#
#  Default command.
#
.DEFAULT:
	@echo "Unknown target $@, try:  make help"

#
#  Local help file.
#
help:
	@echo "   make clean  -  remove all object, module and executable files"
	@echo "   make all    -  make all object, module and executable files"

#
# Dependencies of object files and their compilation rules
#
$(OBJ): %.o: %.f90 
	$(FCC) $(FFLAGS) $(COMPFLAGS) -c $< -o $@ 

#
# Linking stage
#
all: $(OBJ)
	$(FCC) $(FFLAGS) $(OBJ) $(OBJUTIL) -Wl,-rpath=$(MUMPS_ROOT)/lib -Wl,-rpath=$(METIS_ROOT)/lib -Wl,-rpath=$(SCOTCH_ROOT)/lib -o elfe3d 

# clean up
.PHONY: clean

clean:
	rm -f *.o *.mod in/*.vol
