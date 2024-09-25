#!/usr/bin/env python
# coding: utf-8

# ## Configuration of a simulation case (*caseLinear*)
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
from utils import modify_header_file,write_config,write_initial,write_initial_scalar,write_equilibrium,backup_file,restore_file,compile_program,run_program,initialize_variables,read_data_euler,read_data_scalar
import imageio


# ### Setting up the paths
# 
# First, the name of the folder for this test case must be specified:

#Don't forget the bar (/). 
#This directory should have been created prior to the execution to this script, and should also contain an /out folder inside
folder_case="caseLinear/" 


# Then, all the paths are automatically assigned:

#Do not modify the folders and paths below
script_dir = os.path.dirname(os.path.abspath(__file__))
folder_out="out/"
folder_lib="lib"
fname_config="configure.input"
fname_ini="initial.out"
folder_out = os.path.join(script_dir, "../"+folder_case+"/"+folder_out)
folder_case = os.path.join(script_dir, "../"+folder_case)
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
modify_header_file(folder_lib+'/definitions.h', 'NTHREADS', 4)         #number of threads
modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 0)  #System of equations solved
modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)     #Read or not initial data, this should ALWAYS be 1    

# ### Configure the global simulation parameters
# 
# We can set the global simulation parameters as desired:

#Simulation setup
FinalTime = 1.00
DumpTime = 0.05
CFL = 0.4
Order = 7

#Mesh setup
xcells = 100
ycells = 1
zcells = 1
SizeX = 1.0
SizeY = 1.0
SizeZ = 1.0

#Boundary conditions
Face_1 = 1 #-y
Face_2 = 1 #+x
Face_3 = 1 #+y
Face_4 = 1 #-x
Face_5 = 1 #-z
Face_6 = 1 #+z

#Linear transport, only if applicable
u_x = 1.0
u_y = 1.0
u_z = 1.0


# ### Define the initial condition
# 
# To define the initial condition we first need to create the arrays and initialize some variables:

xc, yc, zc, u, uex, *_  = initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ)

for l in range(0,xcells): 
        for m in range(0,ycells): 
            for n in range(0,zcells):
                u[l,m,n]=1.0+np.sin(xc[l,m,n]*2.0*math.pi)                      #imposed as pointwise values for simplicity. For cell averaged data see ordersLinear() in utils.py.
                uex[l,m,n]=1.0+np.sin((xc[l,m,n]-u_x*FinalTime)*2.0*math.pi) 
  
# WRITING CONFIGURATION AND INITIAL DATA

write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, Order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z)      
write_initial_scalar(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u)


# ### Compilation and execution
# 
# The program is compiled and executed:


compile_program()
restore_file(folder_lib+'/definitions.h')
print("Program is running...")
run_program(folder_exe+"./caelum "+folder_case)



# ### Reading data and plotting


files = sorted(glob(folder_out + "/*.out"), key=lambda x: int(re.findall(r'\d+', os.path.basename(x))[0]))
lf = len(files)
lf=len(files)
print(files)

print("Printing figures in folder"+folder_out)

filename = folder_out+"linear_scalar_plot"
fig, ax  = plt.subplots(figsize=(6, 4))
ax.set_xlabel("x") 
ax.set_ylabel("u") 

images = []  # List to hold all the images for the GIF

j=0
u = np.zeros((xcells, ycells, zcells, lf))
for fname in files:
    u = read_data_scalar(u,fname, xcells, ycells, zcells, lf, j)   
    ax.plot(xc[:,0,0],u[:,0,0,j],'o--')
    
    fig2, ax2  = plt.subplots(figsize=(6, 4))
    ax2.plot(xc[:,0,0],u[:,0,0,j],'o-')
    ax2.set_xlabel("x") 
    ax2.set_ylabel("u") 
    image_path = filename + ".png"
    fig2.savefig(image_path,dpi=200)
    images.append(imageio.imread(image_path)) 
    
    j=j+1

ax.plot(xc[:,0,0],u[:,0,0,-1],'o-')
ax.plot(xc[:,0,0],uex[:,0,0],'-k')
ax.set_xlabel("x") 
ax.set_ylabel("u") 

fig.savefig(filename+".png",dpi=300)

gif_path = os.path.join(folder_out, "animation.gif")
imageio.mimsave(gif_path, images, duration=4, loop=0)  # Adjust the duration as needed


    

        