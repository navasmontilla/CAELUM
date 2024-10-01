/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - ibmutils.h

Content:
  -This file contains function prototypes for ibmutils.c

*/


#ifndef IBMUTILS_H
  #define IBMUTILS_H

int assign_image_cells(t_mesh *mesh,t_solid *solids);
int update_ghost_cells(t_sim *sim,t_mesh *mesh,t_solid *solids);
int update_wall_type(t_mesh *mesh,t_solid *solids);

#endif
