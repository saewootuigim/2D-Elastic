ALL: main

FLAGS  = -O2 -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include

OBJECT = global_variable.o my_math.o shape_function.o input.o output.o element_matrix.o find_loading_node.o assemble.o newmark.o main.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: $(OBJECT) chkopts
	-$(FLINKER) $(FLAGS) -o main $(OBJECT) $(PETSC_KSP_LIB)

# post_movie: global_variable.o input.o post_movie.o
# 	ifort $(FFLAGS) $^ -o $@

# post_shannon: global_variable.o input.o post_shannon.o
# 	ifort $(FFLAGS) $^ -o $@
	
# %.o: %.f90
# 	ifort -O3 -c $< -o $@
	
include ${PETSC_DIR}/lib/petsc/conf/test
