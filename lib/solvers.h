/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - solvers.h

Content:
  -This file contains the function prototypes for solvers.c
  
*/

#ifndef SOLVERS_H
  #define SOLVERS_H


  void compute_euler_HLLE(t_wall *wall,double *lambda_max);
  void compute_euler_HLLC(t_wall *wall,double *lambda_max);
  void compute_euler_HLLS(t_wall *wall,double *lambda_max, t_sim *sim);
  void compute_transmissive_euler(t_wall *wall, int wp);
  void compute_solid_euler_hlle(t_wall *wall, double *lambda_max, int wp);
  void compute_euler_Roe(t_wall *wall,double *lambda_max);
  void compute_burgers_flux(t_wall *wall,double *lambda_max);
  void compute_linear_flux(t_wall *wall,double *lambda_max);



#endif

