#!/usr/bin/env python
# coding: utf-8

# ## Configuration of a simulation case (*caseLinear3D*)
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
import pyvista as pv


# ### Setting up the paths
# 
# First, the name of the folder for this test case must be specified:

#Don't forget the bar (/). 
#This directory should have been created prior to the execution to this script, and should also contain an /out folder inside
folder_case="caseLinear3D/" 


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

#Compilation
compile_program()
restore_file(folder_lib+'/definitions.h')

# ### Configure the global simulation parameters
# 
# We can set the global simulation parameters as desired:

#Simulation setup
FinalTime = 1.00
DumpTime = 0.025
CFL = 0.2
Order = 7

#Mesh setup
xcells = 40
ycells = 40
zcells = 40
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
                u[l,m,n]=np.sin(xc[l,m,n]*math.pi)*np.sin(yc[l,m,n]*math.pi)*np.sin(zc[l,m,n]*math.pi)                       #imposed as pointwise values for simplicity. For cell averaged data see ordersLinear() in utils.py.
  
# WRITING CONFIGURATION AND INITIAL DATA

write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, Order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z)      
write_initial_scalar(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u)


# ### Execution
# 
# The program is  executed:

print("Program is running...")
run_program(folder_exe+"./caelum "+folder_case)



# ## RENDERING IN 3D

# Start the virtual framebuffer (Xvfb) to enable off-screen rendering
pv.start_xvfb()

# Find output files
vtk_files = glob(folder_out + "/*.vtk")

# Filter files that contain digits in their filenames
vtk_files_with_digits = [f for f in vtk_files if re.search(r'\d+', os.path.basename(f))]

# Sort the filtered files by the first sequence of digits found in the filenames
files = sorted(vtk_files_with_digits, key=lambda x: int(re.findall(r'\d+', os.path.basename(x))[0]))
lf=len(files)
print(files)

j=0
images = [] 
print("Printing figures in folder"+folder_out)

for fname in files:
    
    vtk_file=fname
    data = pv.read(vtk_file)
    grid = data.compute_cell_sizes().cell_data_to_point_data()
    plotter = pv.Plotter(off_screen=True,window_size=[2200, 1400])

    #opacity = [0.1, 0.1, 0.1, 0.1, 0.9, 0.9, 0.90, 1.0]  
    #plotter.add_volume(grid, scalars="rho", cmap="RdGy_r",opacity=opacity,opacity_unit_distance=0.2,clim=[1.5, 11.5]) #antes 0.1
    
    contour_values = [0.7,0.8,0.9] 
    contours = grid.contour(isosurfaces=contour_values, scalars="U") 
    if contours.n_points > 0:
        plotter.add_mesh(contours, cmap="viridis", scalars="U", opacity=0.7)

    plotter.camera_position = [
        (-1.8, 3.2, 2.1),  # Camera position (x, y, z)
        (0.5, 0.5, 0.5),  # Focal point (center of the object)
        (0, 0, 1),  # View up vector (defines the up direction)
    ]

    
    bounds = [0, 1, 0, 1, 0, 1]  # Adjust to your desired box size
    box = pv.Box(bounds=bounds)
    plotter.add_mesh(box, color="white", line_width=1.0, opacity=0.1) 
    
    plotter.show_axes()

    
    image_path = fname+".png"
    plotter.show(screenshot=image_path)
    img = plotter.screenshot(return_img=True)
    images.append(img)
    
    
    j=j+1
         

gif_path = os.path.join(folder_out, "animation.gif")
imageio.mimsave(gif_path, images, duration=12, loop=0)  # Adjust the duration as needed



    

        