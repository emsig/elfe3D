# choose compiler
#FCC = ifort
FCC = gfortran

# determine operating system
OS = $(shell uname)
# determine architecture
ARCH = $(shell uname -m)
# determine hostname
HOST = $(shell hostname)

# additional directory for shared libraries
LIBDIR = $(HOME)/lib


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


##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
# Section with compiler flags for shared library and compiler options
# Shared libraries:
#   BLAS
#   LAPACK
#   MUMPS
#   SCOTcH
#   METIS

##-------------------------------------------------------------------------------------
##
## for Intel Fortran
##
ifeq (${FCC},ifort)
# for Linux
ifeq (${OS},Linux)
ifeq (${ARCH},x86_64)
# for Intel Fortran
OBJUTIL = -L$(LIBDIR) -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
endif
ifeq (${ARCH},i686)
OBJUTIL = -L/usr/local/lib -L$(LIBDIR)
endif
endif
# for Mac OS X
ifeq (${OS},Darwin)
OBJUTIL = -L/usr/local/lib -L$(LIBDIR) 
endif

# OpenMP library
OMPLIB = 

# for Intel Fortran
FFLAGS = -g -debug all -ftrapuv -traceback -debug inline-debug-info -check all -fpe0 -stand f08 -warn interfaces -warn usage -warn all -diag-enable warn -profile-functions -profile-loops=all -profile-loops-report=2 -assume norealloc_lhs -I/usr/local/include -I$(MKLROOT)/include -qopenmp -mkl=sequential

# compile-time directives (for preprocessor etc.)
# for Linux
ifeq (${OS},Linux)
ifeq (${ARCH},x86_64)
COMPFLAGS = -fpp -DIFORT -DDEBUG
endif
endif
# for Mac OS X
ifeq (${OS},Darwin)
COMPFLAGS = -fpp -DIFORT -DDEBUG
endif
endif

##-------------------------------------------------------------------------------------
##
## for gfortran
##
ifeq (${FCC},gfortran)
# for Linux
ifeq (${OS},Linux)
# for MUMPS
OBJUTIL = -L/usr/local/lib -L$(LIBDIR) -lopenblas -lpthread -lm -lzmumps -lmumps_common -lesmumps -lpord -lmpiseq -L/usr/lib -lscotch -lmetis
endif


# OpenMP library
OMPLIB = -lgomp

FFLAGS = -O2 -march=native -mtune=native -mfma -malign-data=cacheline -std=f2008 -ffree-form -ffree-line-length-none -pedantic -Wall -Wsurprising -Wopenmp-simd -Wvector-operation-performance -ftree-vectorize -ftree-vectorizer-verbose=2 -fopt-info-vec -fopenmp -I/usr/include  -fall-intrinsics

# compile-time directives (for preprocessor etc.)
COMPFLAGS = -x f95-cpp-input -DDEBUG -g
endif
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
# create $(LIBDIR), if it does not exist
#
$(LIBDIR):
	mkdir -p $(LIBDIR)


#
# Dependencies of object files and their compilation rules
#
$(OBJ): %.o: %.f90 
	$(FCC) $(FFLAGS) $(COMPFLAGS) -c $< -o $@ 

#
# Linking stage
#
all: $(OBJ)
	$(FCC) $(FFLAGS) $(OBJ) $(OBJUTIL) -o elfe3d


# clean up
.PHONY: clean
# for Linux
ifeq (${OS},Linux)
clean:
	-\rm *.o *.mod in/*.vol
endif
# for Mac OS X
ifeq (${OS},Darwin)
clean:
	-\rm *.o *.mod in/*.vol
endif

