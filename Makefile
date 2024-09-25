# /*************************************************************************
#
#Authors:
# - Adrián Navas Montilla
# - Isabel Echeverribar

#Copyright (C) 2018-2019 The authors.

#License type: Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Spain (CC BY-NC-ND 3.0 ES https://creativecommons.org/licenses/by-nc-nd/3.0/es/deed.en) under the following terms:

#- Attribution — You must give appropriate credit and provide a link to the license.
#- NonCommercial — You may not use the material for commercial purposes.
#- NoDerivatives — If you remix, transform, or build upon the material, you may not distribute the modified material unless explicit permission of the authors is provided.

#Disclaimer: This software is distributed for research and/or academic purposes, WITHOUT ANY WARRANTY. In no event shall the authors be liable for any claim, damages or other liability, arising from, out of or in connection with the software or the use or other dealings in this software.

#File:
#  - Makefile

#Content:
#  -Make file for set library dependencies and compilation configuration


# Define the compiling core and compilation flags
CC = gcc
DEBUG = 0
OMP = 1
CFLAGS = -Wall

ifeq ($(OMP), 1)
	CFLAGS += -fopenmp
endif

ifeq ($(DEBUG), 1)
	CFLAGS += -g
else
	CFLAGS += -O3
endif

# Define objects and bin file
OBJ = lib/preproc.o lib/ibmutils.o lib/mathutils.o lib/closures.o lib/reconst.o lib/numcore.o lib/solvers.o lib/postproc.o main.o
BIN = caelum

# Rule to construct the exe file
$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm

# Rule to clean
clean:
	$(RM) $(OBJ) $(BIN)

# Rule to generate .o files from .c files
lib/preproc.o: lib/preproc.c
	$(CC) $(CFLAGS) -c -o $@ $<

lib/ibmutils.o: lib/ibmutils.c
	$(CC) $(CFLAGS) -c -o $@ $<

lib/mathutils.o: lib/mathutils.c
	$(CC) $(CFLAGS) -c -o $@ $<

lib/closures.o: lib/closures.c
	$(CC) $(CFLAGS) -c -o $@ $<

lib/reconst.o: lib/reconst.c
	$(CC) $(CFLAGS) -c -o $@ $<

lib/numcore.o: lib/numcore.c
	$(CC) $(CFLAGS) -c -o $@ $<

lib/solvers.o: lib/solvers.c
	$(CC) $(CFLAGS) -c -o $@ $<

lib/postproc.o: lib/postproc.c
	$(CC) $(CFLAGS) -c -o $@ $<

main.o: main.c
	$(CC) $(CFLAGS) -c -o $@ $<