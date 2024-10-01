/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - mathutils.c

Content:
  -This file contains some mathematical functions needed

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>


#include "definitions.h"
#include "structures.h"
#include "mathutils.h"



int vector_product(double *input1, double *input2, double *output){

      output[0]=input1[1]*input2[2]-input1[2]*input2[1];
      output[1]=input1[2]*input2[0]-input1[0]*input2[2];
      output[2]=input1[0]*input2[1]-input1[1]*input2[0];

	return 1;
}

double dot_product(double *input1, double *input2){

      return input1[0]*input2[0]+input1[1]*input2[1]+input1[2]*input2[2];

}

