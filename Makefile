# /*************************************************************************
#
#Authors:
# - Adrián Navas Montilla
# - Isabel Echeverribar
#
#Copyright (C) 2019-2024 The authors.
#
#License type: The 3-Clause BSD License
#
#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
#3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
#This software is provided by the copyright holders and contributors “as is” and any express or implied warranties, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose are disclaimed. In no event shall the copyright holder or contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
#
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