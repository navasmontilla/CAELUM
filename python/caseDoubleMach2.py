#!/usr/bin/env python
# coding: utf-8

# ## Configuration of a simulation case (*caseRP3D*)
# 
# ### Importing required libraries
# 
# Some of these are well known libraries such as *numpy* and *matplotlib* (they can be installed using pip). We also need to import the library *utils* containing predefined functionalities for this software. 

# In[14]:


import os
import math                    
import numpy as np             
import matplotlib.pyplot as plt 
import re
from glob import glob
from utils import modify_header_file,write_config,write_initial,write_equilibrium,backup_file,restore_file,compile_program,run_program,initialize_variables,read_data_euler,write_solid_cells
import imageio
import pyvista as pv


# ### Setting up the paths
# 
# First, the name of the folder for this test case must be specified:

#This test case will run in the folder "caseRP3D/". 
#Don't forget the bar (/). 
#This directory should have been created prior to the execution to this script, and should also contain an /out folder inside
folder_case="caseDoubleMach2/" 


# Then, all the paths are automatically assigned:

#Do not modify the folders and paths below
script_dir = os.path.dirname(os.path.abspath(__file__))
folder_out="out/"
folder_lib="lib"
fname_config="configure.input"
fname_ini="initial.out"
fname_solids="solid_cells.input"
folder_out = os.path.join(script_dir, "../run/"+folder_case+"/"+folder_out)
folder_case = os.path.join(script_dir, "../run/"+folder_case)
folder_lib = os.path.join(script_dir, "../"+folder_lib)
folder_exe = os.path.join(script_dir, "../")

os.makedirs(folder_case, exist_ok=True)
os.makedirs(folder_out, exist_ok=True)

for f in glob(folder_out + "/*.out") + glob(folder_out + "/*.vtk") + glob(folder_out + "/*.png"):
    os.remove(f) 
# ### Setting up the compilation variables (*definitions.h*)
# 
# Here, we can modify those variables that need to be set before compilation and are found in the file *definitions.h*. Don't worry if you mess up things here, a backup of the original file is created before modification and will be restored at the end of this script, after compilation and execution.

#Do not change the line below, it creates a backup of the definitions.h file
backup_file(folder_lib+'/definitions.h')
#Configure the header file for compilation. Add as many lines as desired for the macros you want to modify.
modify_header_file(folder_lib+'/definitions.h', 'NTHREADS', 48)         #number of threads
modify_header_file(folder_lib+'/definitions.h', 'TYPE_REC', 0)
modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 2)  #System of equations solved
modify_header_file(folder_lib+'/definitions.h', 'ST', 0)               #Source term type
modify_header_file(folder_lib+'/definitions.h', 'SOLVER', 0)           #Riemann solver used
modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)     #Read or not initial data, this should ALWAYS be 1
modify_header_file(folder_lib+'/definitions.h', 'WRITE_LIST', 1)  
modify_header_file(folder_lib+'/definitions.h', 'print_POTENTIALTEM', 0)  
modify_header_file(folder_lib+'/definitions.h', 'ALLOW_SOLIDS', 3)  

#Compilation
compile_program()
restore_file(folder_lib+'/definitions.h')     

# ### Configure the global simulation parameters
# 
# We can set the global simulation parameters as desired:

#Simulation setup
FinalTime = 0.2
DumpTime = 0.1
CFL = 0.4 #0.4 in 2D
Order = 7

#Mesh setup
xcells = 300
ycells = 200
zcells = 1
SizeX = 3.0
SizeY = 2.0
SizeZ = SizeY/ycells*zcells

#Boundary conditions
Face_1 = 4 #-y
Face_2 = 3 #+x
Face_3 = 3 #+y
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
sld = np.zeros((xcells, ycells, zcells)) #initialize solid cells


# Among these, we have:
# - the problem variables: ```u, v, w, rho, p, phi```
# - the equilibrium variables (for atmospheric cases): ```ue, ve, we, rhoe, pe```
# 
# Then we can set the initial condition using those variables. To do this, we loop over the three cartesian indexes $(l,m,n)$ and assign the variables, e.g. $\rho(x_l,y_m,z_n)=...$ is set as ```rho[l,m,n]=...```. Cell centers are given by: ```xc[l,m,n]```, ```yc[l,m,n]``` and ```zc[l,m,n]```.

#Initial condition    
for l in range(0,xcells): 
    for m in range(0,ycells):
        for n in range(0,zcells):
            
            if xc[l,m,n] < 0.50:
                p  [l,m,n] = 116.5
                rho[l,m,n] = 8
                u  [l,m,n] = 8.25
                v  [l,m,n] = 0.0
                phi[l,m,n] = 1.0
            else:
                p  [l,m,n] = 1.0
                rho[l,m,n] = 1.4
                u  [l,m,n] = 0.0
                v  [l,m,n] = 0.0
                phi[l,m,n] = 0.0
                

            sld[l,m,n] = 0.5*(-xc[l,m,n] + np.sqrt(3)*yc[l,m,n] + 0.5)
            

# Now, the configuration and initial condition (and equilibrium) files are written: 

write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, Order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z)      
write_initial(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u, v, w, rho, p, phi)
write_solid_cells(folder_case, fname_solids, xcells, ycells, zcells, xc, yc, zc, sld)

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

for fname in files:
    
    u, v, w, rho, p, phi, theta, E = read_data_euler(fname, xcells, ycells, zcells, lf, gamma, j)
                  
    filename = fname+"_doublemach"
    
    xp = xc[:,0,0]     
    yp = yc[0,:,0]      
    X, Y = np.meshgrid(xp, yp)    #matriz de puntos
    Srho=np.transpose(rho[:,:,0,j])
    Spres=np.transpose(p[:,:,0,j])
    Svel=np.transpose(np.sqrt(u[:,:,0,j]**2+w[:,:,0,j]**2))
    Senr=np.transpose(E[:,:,0,j])
    Sphi=np.transpose(phi[:,:,0,j])
    Su=np.transpose(u[:,:,0,j])
    
    lev=[0,1,2,3,4,5,5.5,6,6.5,7,9,9.5,10,11,12,13,14,15,16,17,18,19,20,21,22]
    fig, ax = plt.subplots(figsize=(7, 4))   
    #levels = np.linspace(0, 4, 10)
    #print(levels)
    plot1=ax.contour(X, Y, Srho, levels=lev, colors="k",linewidths=0.2)  
    plot1=ax.contourf(X, Y, Srho, levels=lev, cmap='RdBu')   
    ax.set_xlabel("x") 
    ax.set_ylabel("y") 
    ax.set_xlim([0,3]) 
    ax.set_ylim([0,2]) 
    ax.set_aspect('equal', 'box')
    #plot1.set_clim( 2, 23 )
    # Create colorbar
    cbar = plt.colorbar(plot1)
    cbar.ax.set_title('ρ')
    plt.tight_layout()
    
    # ---- Recta 30° pasando por (0.5,0) ----
    m = np.tan(np.deg2rad(30))
    x_line = np.linspace(0, 3, 500)
    y_line = m * (x_line - 0.5)

    # Dibujar la recta
    ax.plot(x_line, y_line, color='black', linewidth=0.6, zorder=11)

    # Rellenar en blanco por debajo de la recta
    ax.fill_between(x_line, 0, y_line, where=(y_line > 0),
                    color='white', zorder=10)
    
    image_path = filename + ".png"
    fig.savefig(image_path,dpi=250)
    images.append(imageio.imread(image_path)) 
    
    plt.close(fig)
    
    # =====================================================
    # Rotated detail around (0.5, 0) - 30 degrees
    # =====================================================

    x0 = 0.5
    y0 = 0.0
    theta = np.deg2rad(30)

    Lx = 7   # length along rotated direction
    Ly = 3   # width perpendicular

    # Translate
    Xt = X - x0
    Yt = Y - y0

    # Rotate coordinates
    Xr =  Xt*np.cos(theta) + Yt*np.sin(theta)
    Yr = -Xt*np.sin(theta) + Yt*np.cos(theta)

    # Mask rotated rectangle
    mask = (np.abs(Xr) <= Lx/2) & (np.abs(Yr) <= Ly/2)

    Srho_detail = np.where(mask, Srho, np.nan)

    # ---- Create detail figure ----
    fig2, ax2 = plt.subplots(figsize=(6, 4))

    lev=[0,1,2,3,4,5,5.5,6,6.5,7,7.5,8.5,9,9.5,10,10.5,11,12,13,14,15,16,17,18,19,20,21,22]
    plot2 = ax2.contourf(Xr, Yr, Srho_detail, 128, cmap='RdBu')
    ax2.contour(Xr, Yr, Srho_detail, levels=lev, colors="k", linewidths=0.3)

    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_aspect('equal', 'box')

    # Zoom view limits (important!)
    ax2.set_xlim([2.0, 2.8])
    ax2.set_ylim([0, 0.5])

    cbar2 = plt.colorbar(plot2)
    cbar2.ax.set_title('ρ')

    plt.tight_layout()

    # ---- Save detail image with new name ----
    detail_name = fname + "_doublemach_detail"
    detail_path = detail_name + ".png"

    fig2.savefig(detail_path, dpi=250)
    plt.close(fig2)
    
    j=j+1
         

gif_path = os.path.join(folder_out, "animation.gif")
imageio.mimsave(gif_path, images, duration=0.01, loop=0)  # Adjust the duration as needed