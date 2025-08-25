#!/usr/bin/env python
# coding: utf-8

# ## Configuration of a simulation case (*caseCollidingBub*)
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
from utils import modify_header_file,write_config,write_initial,write_equilibrium,backup_file,restore_file,compile_program,run_program,initialize_variables,read_data_euler,write_solid_cells
import imageio

# ### Setting up the paths
# 
# First, the name of the folder for this test case must be specified:

#This test case will run in the folder "caseCollidingBub/". 
#Don't forget the bar (/). 
#This directory should have been created prior to the execution to this script, and should also contain an /out folder inside
folder_case="caseDownwardBubbleSolid/" 


# Then, all the paths are automatically assigned:


#Do not modify the folders and paths below
script_dir = os.path.dirname(os.path.abspath(__file__))
folder_out="out/"
folder_lib="lib"
fname_config="configure.input"
fname_eq="equilibrium.out"
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
modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 2)  #System of equations solved
modify_header_file(folder_lib+'/definitions.h', 'ST', 3)               #Source term type
modify_header_file(folder_lib+'/definitions.h', 'SOLVER', 0)           #Riemann solver used
modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)     #Read or not initial data, this should ALWAYS be 1
modify_header_file(folder_lib+'/definitions.h', 'WRITE_LIST', 1)  
modify_header_file(folder_lib+'/definitions.h', 'ALLOW_SOLIDS', 2)     

#Compilation
compile_program()
restore_file(folder_lib+'/definitions.h')

# ### Configure the global simulation parameters
# 
# We can set the global simulation parameters as desired:

#Simulation setup
FinalTime = 1200.0
DumpTime = 50.0  #for file printing
CFL = 0.45
Order = 7

#Mesh setup
xcells = 200
ycells = 1
zcells = 100
SizeX = 20000.0
SizeY = 1000.0
SizeZ = 10000.0

#Boundary conditions
Face_1 = 4 #-y
Face_2 = 1 #+x
Face_3 = 4 #+y
Face_4 = 1 #-x
Face_5 = 4 #-z
Face_6 = 4 #+z

#Linear transport, only if applicable
u_x = 1.0
u_y = 1.0
u_z = 1.0


# ### Define the initial condition
# 
# To define the initial condition we first need to create the arrays and initialize some variables:

xc, yc, zc, u, v, w, rho, p, phi, ue, ve, we, rhoe, pe = initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ)

sld = np.zeros((xcells, ycells, zcells),dtype=int) #initialize solid cells

# Among these, we have:
# - the problem variables: ```u, v, w, rho, p, phi```
# - the equilibrium variables (for atmospheric cases): ```ue, ve, we, rhoe, pe```
# 
# Then we can set the initial condition using those variables. To do this, we loop over the three cartesian indexes $(l,m,n)$ and assign the variables, e.g. $\rho(x_l,y_m,z_n)=...$ is set as ```rho[l,m,n]=...```. Cell centers are given by: ```xc[l,m,n]```, ```yc[l,m,n]``` and ```zc[l,m,n]```.


#Some auxiliary parameters (must be equal to those in definitions.h)
tt0=300
p0=1.e5
R=287.058
gamma=1.4
g=9.8
rho0=p0/(R*tt0)
aux2=(gamma-1.0)/gamma*g/(R*tt0)
topo = np.zeros((xcells))

#Initial condition and equilibrium state            
for l in range(0,xcells): 
        for m in range(0,ycells): 
            for n in range(0,zcells):
                
                #This is the equilibrium state (in this case, adiabatic equilibrium):
                rhoe[l,m,n]=rho0*(1.0-aux2*zc[l,m,n])**(1.0/(gamma-1.0))
                pe  [l,m,n]=p0*  (1.0-aux2*zc[l,m,n])**(gamma/(gamma-1.0))
                ue  [l,m,n]=0.0
                ve  [l,m,n]=0.0
                we  [l,m,n]=0.0
                
                #This is the initial condition
                xp=10000;
                zp=5000;
                d2=np.sqrt((xc[l,m,n]-xp)*(xc[l,m,n]-xp)+(zc[l,m,n]-zp)*(zc[l,m,n]-zp));

                rc=1000;
                aux1=20.0*min(d2/2.0-rc,0.0)/1000;

                tt=tt0+aux1;
                aux2=(gamma-1.0)/gamma*g/(R*tt0);
                
                rho[l,m,n]=p0/(R*tt)*(1.0-aux2*zc[l,m,n])**(1.0/(gamma-1.0)) #density
                p  [l,m,n]=p0* (1.0-aux2*zc[l,m,n])**(gamma/(gamma-1.0)) #pressure
                u  [l,m,n]=0.0 #x-velocity
                v  [l,m,n]=0.0 #y-velocity
                w  [l,m,n]=0.0 #z-velocity
                phi[l,m,n]=0.0 #solute
                
                topo[l] = (1.0 + np.sin(2*math.pi*xc[l,m,n]/SizeX))*1500
                if  zc[l,m,n] < topo[l]:
                    sld[l, m, n] = 0 #solid cell
                else:
                    sld[l, m, n] = 1 #normal cell


# Now, the configuration and initial condition (and equilibrium) files are written: 

write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, Order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z)
write_equilibrium(folder_case, fname_eq, xcells, ycells, zcells, xc, yc, zc, ue, ve, we, rhoe, pe)       
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

os.remove(folder_out + "list_eq.out")
files = sorted(glob(folder_out + "/*.out"), key=lambda x: int(re.findall(r'\d+', os.path.basename(x))[0]))
#files = glob(folder_out+"/*.out")
lf=len(files)
print(files)

gamma=1.4
j=0
print("Printing figures in folder"+folder_out)
images = []  # List to hold all the images for the GIF
for fname in files:
    
    u, v, w, rho, p, phi, theta, E = read_data_euler(fname, xcells, ycells, zcells, lf, gamma, j)
                  
    filename = fname+"_theta2D"
    
    xp = xc[:,0,0]     
    zp = zc[0,0,:]      
    X, Y = np.meshgrid(xp, zp, indexing="ij")    #matriz de puntos
    
    Sth=theta[:,0,:,j]
    Sth = np.ma.masked_where((Sth < 0) | np.isnan(Sth), Sth)
    
    fig, ax = plt.subplots(figsize=(10, 5))   
    levels = np.linspace(280, 320, 16)
    #print(levels)
    ax.contour(X, Y, Sth, levels=levels,colors="k",linewidths=0.5)  
    plot1=ax.contourf(X, Y, Sth, 200, cmap='RdBu_r', vmin=280, vmax=320)   
    ax.plot(xp, topo, 'k') 
    ax.set_title('Potential temperature')
    ax.set_xlabel("x") 
    ax.set_ylabel("z") 
    ax.set_aspect('equal', 'box')
    #plot1.set_clim( 280, 320 )
    # Create colorbar
    cbar = plt.colorbar(plot1)
    cbar.ax.set_title('θ(K)')
    fig.text(0.15, 0.72, "θ_max="+str(round(np.max(Sth),3)), fontsize=9.5)
    fig.text(0.15, 0.68, "θ_min="+str(round(np.min(Sth),3)), fontsize=9.5)
    
    image_path = filename + ".png"
    fig.savefig(image_path,dpi=400)
    images.append(imageio.imread(image_path)) 
    
    j=j+1
         
gif_path = os.path.join(folder_out, "animation.gif")
imageio.mimsave(gif_path, images, duration=8, loop=0)  # Adjust the duration as needed
