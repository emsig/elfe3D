# choose compiler
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
# Section with compiler flags for shared library and compiler options
##
## for gfortran
##
ifeq (${FCC},gfortran)
ifeq (${OS},Linux)
# link MUMPS libaries here!
OBJUTIL = -L/usr/local/lib -L$(LIBDIR) -lpthread -lm -lzmumps -lmumps_common -lesmumps -lpord -lmpiseq -L/usr/lib -lscotch -lmetis -lopenblas
endif

# OpenMP library
OMPLIB = -lgomp

# compiler flags
FFLAGS = -O2 -march=native -mtune=native -mfma -malign-data=cacheline -std=f2008 -ffree-form -ffree-line-length-none -ftree-vectorize -ftree-vectorizer-verbose=2 -fopt-info-vec -fopenmp -I/usr/include -fall-intrinsics
# for warnings you can include here: 
# -pedantic -Wall -Wopenmp-simd -Wvector-operation-performance

# compile-time directives (for preprocessor etc.)
COMPFLAGS = -x f95-cpp-input
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

ifeq (${OS},Linux)
clean:
	-\rm *.o *.mod in/*.vol
endif
