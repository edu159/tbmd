IDIR=include
CC=g++ -O2 #-Wall -pedantic
CFLAGS=

ODIR=build
SRCDIR=src
BINDIR=bin

LIBS= -llapack -lgsl -lblas


_SRC = main.cpp File.cpp OutputDataFile.cpp ConfigFile.cpp InputDataFile.cpp MovieFile.cpp AttractiveEnergy.cpp  common.cpp  Energy.cpp Hamiltonian.ccp  MDSimulation.cpp  Particle.cpp  ParticleSystem.cpp  RepulsiveEnergy.cpp
SRC = $(patsubst %,$(SRCDIR)/%,$(_SRC))

_DEPS = File.h OutputDataFile.h ConfigFile.h InputDataFile.h MovieFile.h AttractiveEnergy.h  common.h  Energy.h  Hamiltonian.h  MDSimulation.h  Particle.h  ParticleSystem.h  RepulsiveEnergy.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o File.o OutputDataFile.o ConfigFile.o InputDataFile.o MovieFile.o AttractiveEnergy.o  common.o  Energy.o  Hamiltonian.o  MDSimulation.o  Particle.o  ParticleSystem.o  RepulsiveEnergy.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

default: mkdirs tbmd

mkdirs:
	@mkdir -p $(ODIR)
	@mkdir -p $(SRCDIR)
	@mkdir -p $(BINDIR)

tbmd: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) 
	@mv tbmd $(BINDIR)

$(OBJ): $(ODIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -I $(IDIR) -c $< -o $@


.PHONY: clean dist

clean:
	@rm -rf $(ODIR) $(BINDIR) $(ODIR)

dist:
	tar -zcvf TBMD.tar.gz * --verbose

