########################################
# Replace the main_program below with the name of
# your main program without the ".f90" extension
# replace the module_name below with your module
# name without the ".f90" extension

# PROGRAM NAME (no .f90)
main = main

# MODULE NAME (no .f90)
mod0 = module_values
mod1 = module_ode
mod2 = module_interpolate
########################################
cmplr = gfortran
flags = -g -Wall

objects1 = $(mod0).o $(mod1).o	$(mod2).o		$(main).o

$(main)		:	$(objects1)
				$(cmplr) $(flags) -o $(main) $(objects1) -llapack -lblas
# -llapack -lblas flags to access outside libraries

$(mod0).o : $(mod0).f90
				$(cmplr) -c $(mod0).f90

$(mod1).o	:	$(mod0).f90		$(mod1).f90
				$(cmplr) -c $(mod1).f90

$(mod2).o :	$(mod0).f90		$(mod2).f90
				$(cmplr) -c $(mod2).f90

$(main).o	:	$(main).f90		$(mod0).f90		$(mod1).f90		$(mod2).f90
				$(cmplr) -c $(main).f90

clean			:
				rm -f *.mod *.pcl *.pc *.o $(main) *.inc *.vo *.d *.out
