###################################################################


FC= ifort 


#.. Set executable name
PROGRAM= efimov.x

###################################################################


#.. Set compiler flags
FLAGS= -O3

LIBRARIES = -lblas -llapack 


############################################################

OBJECTS = constants.o gauleg.o universal_knot.o efimov_sj.o efimovham_sj2.o twobody_potential_sj.o model_potential_sj.o bder.o bget.o  bsplvb.o setup_sj.o

# Linking to executable

$(PROGRAM): $(OBJECTS)
	$(FC)  $(LIBRARIES) $(FLAGS)  -o $(PROGRAM) $(OBJECTS)

# Module

constants.mod: constants.o constants_sj.f90
	$(FC) $(FLAGS) -c constants_sj.f90


nrtype.mod: nrtype.o nrtype.f90
	$(FC) $(FLAGS) -c nrtype.f90
nrtype.o: nrtype.f90
	$(FC) $(FLAGS) -c nrtype.f90

nr.mod: nr.o nr.f90
	$(FC) $(FLAGS) -c nr.f90 
nr.o: nrtype.o nr.f90
	$(FC) $(FLAGS) -c nr.f90
nrutil.mod: nrutil.o nrutil.f90
	$(FC) $(FLAGS) -c nrutil.f90
nrutil.o: nrtype.o nrutil.f90
	$(FC) $(FLAGS) -c nrutil.f90

# Object files

%.o: %.f90
	$(FC) $(FLAGS) -c $<

# remove files

clean:
	rm $(OBJECTS)  $(PROGRAM) constants.mod
# End of the makefile


