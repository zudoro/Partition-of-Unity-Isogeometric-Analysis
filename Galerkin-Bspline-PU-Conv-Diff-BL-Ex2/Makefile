#########################################################

# REPLACE the main_prg_name below with the name of
# your main program without the ".f90" extension
# REPLACE the module_name below with your module 
# name without the ".f90" extension

# PROGRAM NAME (no .f90)
main = MAIN
exec = test_run

# MODULE NAME  (no .f90)
mod1=NEWTYPE
mod2=GLBVAR
mod3=KNOT_HANDLING
mod4=GSQUAD
mod5=LUDECOMPOSITION
mod6=NURBS
mod7=PATCH_MAPPING
mod8=GEOMETRY
mod9=PLOT
mod10=NURBS_BASIS
mod11=LOADFUNCTION
mod12=INTEGRATION
mod13=BOUNDARY
mod14=ERRORESTIMATE
#########################################################

cmplr = gfortran
flag = -ffixed-line-length-none -ffree-line-length-none

# objects1 = $(mod1).o $(mod2).o $(mod3).o $(mod4).o $(mod5).o $(mod6).o $(mod7).o $(mod8).o $(mod9).o $(mod10).o $(mod11).o $(mod12).o $(mod13).o $(mod14).o $(main).o
# objects2 = $(mod1).f90 $(mod2).f90 $(mod3).f90 $(mod4).f90 $(mod5).f90 $(mod6).f90 $(mod7).f90 $(mod8).f90 $(mod9).f90 $(mod10).f90 $(mod11).f90 $(mod12).f90 $(mod13).f90 $(mod14).f90
# 
# $(main)   : $(objects1)
# 	$(cmplr) $(flag) -o $(exec) $(objects1)
# 
# $(mod1).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod2).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod3).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod4).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod5).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod6).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod7).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod8).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod9).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod10).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod11).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod12).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod13).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(mod14).o    : $(objects2)
# 	$(cmplr) $(flag) -c $(objects2)
# 
# $(main).o  : $(main).f90 $(objects2)
# 	$(cmplr) $(flag) -c $(objects2) $(main).f90
# 
# clean :
# 	rm -f *.mod *.o
# 
# cleanall :
# 	rm -f *.mod *.o *~ $(exec)

objects1 = $(mod1).o $(mod2).o $(mod3).o $(mod4).o $(mod5).o $(mod6).o $(mod7).o $(mod8).o $(mod9).o $(mod10).o $(mod11).o $(mod12).o $(mod13).o $(mod14).o $(main).o
objects2 = $(mod1).f90 $(mod2).f90 $(mod3).f90 $(mod4).f90 $(mod5).f90 $(mod6).f90 $(mod7).f90 $(mod8).f90 $(mod9).f90 $(mod10).f90 $(mod11).f90 $(mod12).f90 $(mod13).f90 $(mod14).f90

$(main)   : $(objects1)
	$(cmplr) $(flag) -o $(exec) $(objects1)

$(mod1).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod2).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod3).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod4).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod5).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod6).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod7).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod8).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod9).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod10).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod11).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod12).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod13).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(mod14).o    : $(objects2)
	$(cmplr) $(flag) -c $(objects2) -O3

$(main).o  : $(main).f90 $(objects2)
	$(cmplr) $(flag) -c $(objects2) $(main).f90 -O3

clean :
	rm -f *.mod *.o

cleanall :
	rm -f *.mod *.o *~ $(exec)