# Makefile for source/Simulation/SimulationMain/Cellular

Simulation += Simulation_data.o \
	sim_ranmar.o User_interface.o User_initBlobCell.o User_updateEnvelopeComposition.o \
	User_calculateAtmosphere.o
Simulation_init.o : Simulation_data.o
Simulation_initBlock.o : Simulation_data.o
User_initBlobCell.o : Simulation_data.o sim_ranmar.o

