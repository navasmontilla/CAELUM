#!/usr/bin/env python
# coding: utf-8

# ## Configuration of a simulation case (*caseRM3D*)
# 
# ### Importing required libraries
# 
# Some of these are well known libraries such as *numpy* and *matplotlib* (they can be installed using pip). We also need to import the library *utils* containing predefined functionalities for this software. 

import os
import math                    
import numpy as np             
import matplotlib.pyplot as plt 
import re
from glob import glob
from utils import modify_header_file,write_config,write_initial,write_equilibrium,backup_file,restore_file,compile_program,run_program,initialize_variables,read_data_euler
import imageio


# ### Setting up the paths
# 
# First, the name of the folder for this test case must be specified:

#This test case will run in the folder "caseShockBub/". 
#Don't forget the bar (/). 
#This directory should have been created prior to the execution to this script, and should also contain an /out folder inside
folder_case="caseRM3D/" 


# Then, all the paths are automatically assigned:

#Do not modify the folders and paths below
script_dir = os.path.dirname(os.path.abspath(__file__))
folder_out="out/"
folder_lib="lib"
fname_config="configure.input"
fname_ini="initial.out"
folder_out = os.path.join(script_dir, "../run/"+folder_case+"/"+folder_out)
folder_case = os.path.join(script_dir, "../run/"+folder_case)
folder_lib = os.path.join(script_dir, "../"+folder_lib)
folder_exe = os.path.join(script_dir, "../")

os.makedirs(folder_case, exist_ok=True)
os.makedirs(folder_out, exist_ok=True)
# ### Setting up the compilation variables (*definitions.h*)
# 
# Here, we can modify those variables that need to be set before compilation and are found in the file *definitions.h*. Don't worry if you mess up things here, a backup of the original file is created before modification and will be restored at the end of this script, after compilation and execution.

#Do not change the line below, it creates a backup of the definitions.h file
backup_file(folder_lib+'/definitions.h')
#Configure the header file for compilation. Add as many lines as desired for the macros you want to modify.
modify_header_file(folder_lib+'/definitions.h', 'NTHREADS', 60)         #number of threads
modify_header_file(folder_lib+'/definitions.h', 'TYPE_REC', 0)
modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 2)  #System of equations solved
modify_header_file(folder_lib+'/definitions.h', 'ST', 0)               #Source term type
modify_header_file(folder_lib+'/definitions.h', 'SOLVER', 1)           #Riemann solver used
modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)     #Read or not initial data, this should ALWAYS be 1   
modify_header_file(folder_lib+'/definitions.h', 'WRITE_LIST', 0)       
modify_header_file(folder_lib+'/definitions.h', 'print_VELOCITY', 0)     
modify_header_file(folder_lib+'/definitions.h', 'print_POTENTIALTEM', 0)     
modify_header_file(folder_lib+'/definitions.h', 'print_PRESSURE', 1)   
modify_header_file(folder_lib+'/definitions.h', 'print_RHO', 1)    

#Compilation
compile_program()
restore_file(folder_lib+'/definitions.h')

# ### Configure the global simulation parameters
# 
# We can set the global simulation parameters as desired:

#Simulation setup
FinalTime = 5.0
DumpTime = 0.2  #for file printing
CFL = 0.2
Order = 7

#Mesh setup
xcells = 20
ycells = 20
zcells = 20
SizeX = 1.0
SizeY = 1.0
SizeZ = 1.0

#Boundary conditions
Face_1 = 1 #-y
Face_2 = 3 #+x
Face_3 = 1 #+y
Face_4 = 3 #-x
Face_5 = 1 #-z
Face_6 = 1 #+z

#Linear transport, only if applicable
u_x = 1.0
u_y = 1.0
u_z = 1.0


# ### Define the initial condition
# 
# To define the initial condition we first need to create the arrays and initialize some variables:

xc, yc, zc, u, v, w, rho, p, phi, ue, ve, we, rhoe, pe = initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ)


# Among these, we have:
# - the problem variables: ```u, v, w, rho, p, phi```
# - the equilibrium variables (for atmospheric cases): ```ue, ve, we, rhoe, pe```
# 
# Then we can set the initial condition using those variables. To do this, we loop over the three cartesian indexes $(l,m,n)$ and assign the variables, e.g. $\rho(x_l,y_m,z_n)=...$ is set as ```rho[l,m,n]=...```. Cell centers are given by: ```xc[l,m,n]```, ```yc[l,m,n]``` and ```zc[l,m,n]```.

#Initial condition and equilibrium state            
for l in range(0,xcells): 
        for m in range(0,ycells): 
            for n in range(0,zcells):
                u[l,m,n] = 0.0
                v[l,m,n] = 0.0
                w[l,m,n] = 0.0

                # Compute conditions
                if xc[l,m,n] < 0.47:
                    p[l,m,n] = 1.86  # 10.956
                    rho[l,m,n] = 1.93  # 3.936
                    u[l,m,n] = 0.26  # 2.7245
                    phi[l,m,n] = 0.0
                elif xc[l,m,n] < 0.59 + 0.04 * np.cos(2.0 * math.pi * yc[l,m,n]) + 0.04 * np.cos(2.0 * math.pi * zc[l,m,n]):
                    p[l,m,n] = 0.77  # 1.0
                    rho[l,m,n] = 1.0
                    u[l,m,n] = -0.47  # 0.0
                    v[l,m,n] = 0.0
                    phi[l,m,n] = 0.0
                else:
                    p[l,m,n] = 0.77  # 1.0
                    rho[l,m,n] = 4.88  # 0.50  # 0.1
                    u[l,m,n] = -0.47  # 0.0
                    phi[l,m,n] = 1.0


# Now, the configuration and initial condition (and equilibrium) files are written: 

write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, Order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z)      
write_initial(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u, v, w, rho, p, phi)


# ### Execution
# 
# The program is  executed:

print("Program is running...")
run_program(folder_exe+"./caelum "+folder_case)



# ### Reading data and plotting
# 
# To read data, we use the function ```read_data_euler()``` which gives as output the numerical solution for all time levels, e.g.  $\rho(x_l,y_m,z_n,t_j)$ is ```rho[l,m,n,j]```, where ```j``` is the time level.
# 
# This can be customized for each particular case. 

files = sorted(glob(folder_out + "/*.out"), key=lambda x: int(re.findall(r'\d+', os.path.basename(x))[0]))
lf = len(files)
lf=len(files)
print(files)

gamma=1.4
j=0
print("Printing figures in folder"+folder_out)

images = []  # List to hold all the images for the GIF

# for fname in files:
    
    # u, v, w, rho, p, phi, theta, E = read_data_euler(fname, xcells, ycells, zcells, lf, gamma, j)
                  
    # filename = fname+"_shockbubble"
    
    # xp = xc[:,0,0]     
    # yp = yc[0,:,0]      
    # X, Y = np.meshgrid(xp, yp)    #matriz de puntos
    # Srho=np.transpose(rho[:,:,0,j])

    
    # fig, ax = plt.subplots(figsize=(10, 5))   
    # levels = np.linspace(0, 4, 10)
    # plot1=ax.contour(X, Y, Srho, levels=levels,colors="k",linewidths=0.2)  
    # plot1=ax.contourf(X, Y, Srho, cmap='plasma')   
    # ax.set_title('Density')
    # ax.set_xlabel("x") 
    # ax.set_ylabel("y") 
    # ax.set_aspect('equal', 'box')
    # plot1.set_clim( 0, 4 )
    # Create colorbar
    # cbar = plt.colorbar(plot1)
    # cbar.ax.set_title('ρ')
    # fig.text(0.15, 0.72, "ρ_max="+str(round(np.max(Srho),3)), fontsize=9.5)
    # fig.text(0.15, 0.68, "ρ_min="+str(round(np.min(Srho),3)), fontsize=9.5)
    
    # image_path = filename + ".png"
    #fig.savefig(image_path,dpi=400)
    #images.append(imageio.imread(image_path)) 
    
    # j=j+1
         

#gif_path = os.path.join(folder_out, "animation.gif")
#imageio.mimsave(gif_path, images, duration=8, loop=0)  # Adjust the duration as needed