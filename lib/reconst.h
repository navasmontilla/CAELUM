/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - reconst.h

Content:
  -This file contains the function prototypes for reconst.c
  
*/



#ifndef RECONST_H
  #define RECONST_H


  double weno3L(double *phi);
  double weno3R(double *phi);
  double weno5L(double *phi);
  double weno5R(double *phi);
  double weno7L(double *phi);
  double weno7R(double *phi);



#endif

