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
  - declaration.h

Content:
  -This header file contains the declaration of all the functions which are implemented in the code
  They are grouped by libraries in this file, and distributed in different C files in the lib folder.

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

