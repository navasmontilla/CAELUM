#!/usr/bin/env python
# coding: utf-8

# ## Configuration of a simulation case (*caseShockBub*)
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

#Don't forget the bar (/). 
#This directory should have been created prior to the execution to this script, and should also contain an /out folder inside
folder_case="caseRP/" 


# Then, all the paths are automatically assigned:

#Do not modify the folders and paths below
script_dir = os.path.dirname(os.path.abspath(__file__))
folder_out="out/"
folder_lib="lib"
folder_exact="autotest/caseRPs/exact/"
fname_config="configure.input"
fname_ini="initial.out"
folder_out = os.path.join(script_dir, "../"+folder_case+"/"+folder_out)
folder_case = os.path.join(script_dir, "../"+folder_case)
folder_lib = os.path.join(script_dir, "../"+folder_lib)
folder_exact = os.path.join(script_dir, "../"+folder_exact)
folder_exe = os.path.join(script_dir, "../")

os.makedirs(folder_case, exist_ok=True)
os.makedirs(folder_out, exist_ok=True)

for f in glob(folder_out + "/*.out") + glob(folder_out + "/*.vtk") + glob(folder_out + "/*.png"):
    os.remove(f) 
# ### Setting up the compilation variables (*definitions.h*)
# 
# Here, we can modify those variables that need to be set before compilation and are found in the file *definitions.h*. Don't worry if you mess up things here, a backup of the original file is created before modification and will be restored at the end of this script, after compilation and execution.

# In[17]:


#Do not change the line below, it creates a backup of the definitions.h file
backup_file(folder_lib+'/definitions.h')
#Configure the header file for compilation. Add as many lines as desired for the macros you want to modify.
modify_header_file(folder_lib+'/definitions.h', 'NTHREADS', 4)         #number of threads
modify_header_file(folder_lib+'/definitions.h', 'TYPE_REC', 0)
modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 2)  #System of equations solved
modify_header_file(folder_lib+'/definitions.h', 'ST', 0)               #Source term type
modify_header_file(folder_lib+'/definitions.h', 'SOLVER', 0)           #Riemann solver used
modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)     #Read or not initial data, this should ALWAYS be 1    

# ### Configure the global simulation parameters
# 
# We can set the global simulation parameters as desired:

# In[18]:


#Simulation setup
DumpTime = 0.5
CFL = 0.4
Order = 5

#Mesh setup
xcells = 200
ycells = 1
zcells = 1
SizeX = 1.0
SizeY = 1.0
SizeZ = 1.0

#Boundary conditions
Face_1 = 3 #-y
Face_2 = 3 #+x
Face_3 = 3 #+y
Face_4 = 3 #-x
Face_5 = 3 #-z
Face_6 = 3 #+z

#Linear transport, only if applicable
u_x = 1.0
u_y = 1.0
u_z = 1.0


# ### Define the initial condition
# 
# To define the initial condition we first need to create the arrays and initialize some variables:

# In[19]:


xc, yc, zc, u, v, w, rho, p, phi, ue, ve, we, rhoe, pe = initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ)


# Among these, we have:
# - the problem variables: ```u, v, w, rho, p, phi```
# - the equilibrium variables (for atmospheric cases): ```ue, ve, we, rhoe, pe```
# 
# Then we can set the initial condition using those variables. To do this, we loop over the three cartesian indexes $(l,m,n)$ and assign the variables, e.g. $\rho(x_l,y_m,z_n)=...$ is set as ```rho[l,m,n]=...```. Cell centers are given by: ```xc[l,m,n]```, ```yc[l,m,n]``` and ```zc[l,m,n]```.

# In[20]:

case = 2

if case==1: #RP1 Steady
    FinalTime = 0.011
    for l in range(0,xcells): 
        if (xc[l,:,:]<0.5):
            rho[l,:,:]=1.0
            p  [l,:,:]=1.0
            u  [l,:,:]=0.0
            phi[l,:,:]=1.0
        else:
            rho[l,:,:]=1.0
            p  [l,:,:]=1.0
            u  [l,:,:]=0.0
            phi[l,:,:]=0.0 


if case==2: #RP2 Sod shock
    FinalTime = 0.2
    DumpTime = 0.01
    for l in range(0,xcells): 
        if (xc[l,:,:]<0.5):
            rho[l,:,:]=1.0
            p  [l,:,:]=1.0
            u  [l,:,:]=0.0
            phi[l,:,:]=1/(1.6-1) 
        else:
            rho[l,:,:]=0.125
            p  [l,:,:]=0.1
            u  [l,:,:]=0.0
            phi[l,:,:]=1/(1.2-1)
    exactS  = np.loadtxt(folder_exact+"RP1.txt") 
    
if case==3: #RP3
    FinalTime = 0.035
    DumpTime = 0.002
    for l in range(0,xcells): 
        if (xc[l,:,:]<0.4):
            rho[l,:,:]=5.99924
            p  [l,:,:]=460.894
            u  [l,:,:]=19.5975
            phi[l,:,:]=1.0
        else:
            rho[l,:,:]=5.99242
            p  [l,:,:]=46.0950
            u  [l,:,:]=-6.19633
            phi[l,:,:]=0.0
    exactS  = np.loadtxt(folder_exact+"RP3.txt")  
    
if case==4: #RP4
    modify_header_file(folder_lib+'/definitions.h', 'MULTICOMPONENT', 1)
    FinalTime = 0.16
    DumpTime = 0.01
    for l in range(0,xcells): 
        if (xc[l,:,:]<0.5):
            rho[l,:,:]=1.0
            p  [l,:,:]=1.0
            u  [l,:,:]=0.0
            phi[l,:,:]=1/(1.4-1) #1/(gamma-1) for gamma formulation
        else:
            rho[l,:,:]=0.125
            p  [l,:,:]=0.1
            u  [l,:,:]=0.0
            phi[l,:,:]=1/(1.6-1)
    exactSm = np.loadtxt(folder_exact+"RP1_multi.txt")   
    
    
if case==5: #RP5 
    FinalTime = 0.011
    DumpTime = 0.001
    for l in range(0,xcells): 
        if (xc[l,:,:]<0.5):
            rho[l,:,:]=1.0
            p  [l,:,:]=1000.0
            u  [l,:,:]=0.0
            phi[l,:,:]=1.0
        else:
            rho[l,:,:]=1.0
            p  [l,:,:]=0.01
            u  [l,:,:]=0.0
            phi[l,:,:]=0.0
    exactS  = np.loadtxt(folder_exact+"RP2.txt")  


# Now, the configuration and initial condition (and equilibrium) files are written: 

# In[21]:


write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, Order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z)      
write_initial(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u, v, w, rho, p, phi)


# ### Compilation and execution
# 
# The program is compiled and executed:

# In[22]:


compile_program()
print("Program is running...")
run_program(folder_exe+"./exehow3d "+folder_case)
restore_file(folder_lib+'/definitions.h')


# ### Reading data and plotting
# 
# To read data, we use the function ```read_data_euler()``` which gives as output the numerical solution for all time levels, e.g.  $\rho(x_l,y_m,z_n,t_j)$ is ```rho[l,m,n,j]```, where ```j``` is the time level.
# 
# This can be customized for each particular case. 

# In[23]:


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
                  
    filename = fname+"_RP"+str(case)
    
    fig, ax  = plt.subplots(2,2,figsize=(14, 8))
    ax[0,0].plot(xc[:,0,0],rho[:,0,0,j],'o-')
    if 'exactS' in locals():
        ax[0,0].plot(exactS[:,0],exactS[:,1],'k-') 
    if 'exactSm' in locals():
        ax[0,0].plot(exactSm[:,0],exactSm[:,1],'k-') 
    ax[0,0].set_xlabel("x") 
    ax[0,0].set_ylabel("density") 

    ax[0,1].plot(xc[:,0,0],p[:,0,0,j],'o-')
    if 'exactS' in locals():
        ax[0,1].plot(exactS[:,0],exactS[:,2],'k-') 
    ax[0,1].set_xlabel("x") 
    ax[0,1].set_ylabel("pressure") 

    ax[1,0].plot(xc[:,0,0],u[:,0,0,j],'o-') 
    if 'exactS' in locals():
        ax[1,0].plot(exactS[:,0],exactS[:,3],'k-') 
    ax[1,0].set_xlabel("x") 
    ax[1,0].set_ylabel("velocity u") 

    ax[1,1].plot(xc[:,0,0],phi[:,0,0,j],'o-') 
    ax[1,1].set_xlabel("x") 
    ax[1,1].set_ylabel("phi")  
    
    image_path = filename + ".png"
    fig.savefig(image_path,dpi=200)
    images.append(imageio.imread(image_path)) 
    
    j=j+1
         

gif_path = os.path.join(folder_out, "animation.gif")
imageio.mimsave(gif_path, images, duration=8, loop=0)  # Adjust the duration as needed