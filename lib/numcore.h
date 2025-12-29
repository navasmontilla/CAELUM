/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - numcore.h

Content:
  -This file contains all the function prototypes for numcore.c
  
*/


#ifndef NUMCORE_H
  #define NUMCORE_H

  void update_cell(t_mesh *mesh, t_sim *sim);
  void update_cellK1(t_mesh *mesh, t_sim *sim);
  void update_cellK2(t_mesh *mesh, t_sim *sim);
  void update_cellK3(t_mesh *mesh, t_sim *sim);

  int equilibrium_reconstruction(t_mesh *mesh, t_sim *sim);
  int compute_fluxes(t_mesh *mesh, t_sim *sim);
  void compute_transport(t_wall *wall);
  void compute_source(t_mesh *mesh);

  int update_cell_boundaries(t_mesh *mesh);
  int update_dt(t_mesh *mesh,t_sim *sim);

  void mass_calculation(t_mesh *mesh, t_sim *sim);
  void energy_calculation(t_mesh *mesh, t_sim *sim);
  void tke_calculation(t_mesh *mesh, t_sim *sim);
  
  void update_solution(t_mesh *mesh, t_sim *sim, t_solid *solids, int rk_steps);
  
  void positivity_fix(double nvar, t_cell *cell, double *UL, double *UR, double zL, double zR);
  void positivity_fix_loop(t_mesh *mesh, t_sim *sim);


#endif

