/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - closures.h

Content:
  -This code file contains the function prototypes for closures.c

*/


#ifndef CLOSURES_H
  #define CLOSURES_H

double energy_from_pressure(double gm, double p, double u, double v, double w, double rho, double z);
double pressure_from_energy(double gm, double E, double u, double v, double w, double rho, double z);

#endif
