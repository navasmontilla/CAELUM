/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2018-2019 The authors.

License type: Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Spain (CC BY-NC-ND 3.0 ES https://creativecommons.org/licenses/by-nc-nd/3.0/es/deed.en) under the following terms:

- Attribution — You must give appropriate credit and provide a link to the license.
- NonCommercial — You may not use the material for commercial purposes.
- NoDerivatives — If you remix, transform, or build upon the material, you may not distribute the modified material unless explicit permission of the authors is provided.

Disclaimer: This software is distributed for research and/or academic purposes, WITHOUT ANY WARRANTY. In no event shall the authors be liable for any claim, damages or other liability, arising from, out of or in connection with the software or the use or other dealings in this software.

File:
  - mathutils.c

Content:
  -This code file contains all the auxiliar mathematical functions needed

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

