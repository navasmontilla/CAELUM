/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - postproc.h

Content:
  -This file contains function prototypes for postproc.c
  
*/


#ifndef POSTPROC_H
  #define POSTPROC_H

  int write_vtk(t_mesh *mesh,char *filename);
  int write_list(t_mesh *mesh,char *filename);
  int write_list_eq(t_mesh *mesh,char *filename);
  int write_geo_vtk(t_mesh *mesh, char *filename);
  void screen_info(t_mesh *mesh, t_sim *sim);

#endif

