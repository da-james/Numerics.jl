########################################
# Replace the main_program below with the name of
# your main program without the ".f90" extension
# replace the module_name below with your module
# name without the ".f90" extension

# PROGRAM NAME (no .f90)
main = main

# MODULE NAME (no .f90)
mod1 = module_ode
########################################
cmplr = gfortran
flags = -g -Wall

objects1 = $(mod1).o	$(main).o

$(main)		:	$(objects1)
				$(cmplr) $(flags) -o $(main) $(objects1) -llapack -lblas
# -llapack -lblas flags to access outside libraries

$(mod1).o	:	$(mod1).f90
				$(cmplr) -c $(mod1).f90

$(main).o	:	$(main).f90		$(mod1).f90
				$(cmplr) -c $(main).f90

clean			:
				rm -f *.mod *.pcl *.pc *.o $(main) *.inc *.vo *.d *.out
