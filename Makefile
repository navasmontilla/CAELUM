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


# Definir el compilador y las banderas básicas de compilación
CC = gcc
DEBUG = 0
CFLAGS =  -Wall  -fopenmp

ifeq ($(DEBUG), 1)
	CFLAGS += -g
else
	CFLAGS += -O3
endif

# Definir los archivos objeto y el archivo binario
OBJ = lib/preproc.o lib/ibmutils.o lib/mathutils.o lib/closures.o lib/reconst.o lib/numcore.o lib/solvers.o lib/postproc.o ehow3d.o
BIN = exehow3d

# Regla para construir el ejecutable
$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm

# Regla para limpiar los archivos generados
clean:
	$(RM) $(OBJ) $(BIN)

# Reglas individuales para compilar cada archivo .c en su correspondiente archivo .o
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

ehow3d.o: ehow3d.c
	$(CC) $(CFLAGS) -c -o $@ $<