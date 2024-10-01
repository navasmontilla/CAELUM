/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - preproc.h

Content:
  -This file contains the function prototypes for preproc.c

*/



#ifndef PREPROC_H
  #define PREPROC_H

  int create_mesh(t_mesh *mesh,t_sim *sim);
  int update_initial(t_mesh *mesh, t_sim *sim, const char *folder_path);
  int read_initial(t_mesh *mesh, t_sim *sim, const char *folder_path);
  int assign_wall_type(t_mesh *mesh);
  int assign_cell_type(t_mesh *mesh,t_solid *solids);
  int update_stencils(t_mesh *mesh,t_sim *sim);
  int read_solids(t_mesh *mesh,t_solid *solids, const char *folder_path);
  int read_solids_txt(t_mesh *mesh,t_solid *solids, const char *folder_path);
  void set_velocity_field(t_mesh *mesh, t_sim *sim);  //ESTO LO QUITAREMOS, JUNTO CON LA MACRO "LINEAR_TRANSPORT" y sus ifs asociados
  void set_velocity(t_mesh *mesh, t_sim *sim);
  void read_config(t_mesh *mesh, t_sim *sim, const char *folder_path);
  void print_info(t_mesh *mesh, t_sim *sim, const char *folder_path);

#endif
