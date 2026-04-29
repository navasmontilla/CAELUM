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
folder_case="caseDownwardBubbleIBM/" 


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
modify_header_file(folder_lib+'/definitions.h', 'ST', 2)               #Source term type
modify_header_file(folder_lib+'/definitions.h', 'SOLVER', 0)           #Riemann solver used
modify_header_file(folder_lib+'/definitions.h', 'TYPE_REC', 2)
modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)     #Read or not initial data, this should ALWAYS be 1
modify_header_file(folder_lib+'/definitions.h', 'WRITE_LIST', 1)  
modify_header_file(folder_lib+'/definitions.h', 'ALLOW_SOLIDS', 3)     

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
xcells = 400
ycells = 1
zcells = 200
SizeX = 20000.0
SizeY = SizeX/xcells
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

sld = np.zeros((xcells, ycells, zcells)) #initialize solid cells

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
               
                #sld[l,m,n] =  zc[l,m,n]-3000  
                sld[l,m,n] = np.sqrt((xc[l,m,n]-xp)**2 + (zc[l,m,n]+2500.0)**2) - 5000                
                #sld[l,m,n] =  1.0/np.sqrt(1+0.1**2)*(-xc[l,m,n] + 10*zc[l,m,n] - 5000) 


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

xc_circ = 10000.0
zc_circ = -2500.0
R_circ  = 5000.0

def sdf_circle(x, z, xc, zc, R):
    return np.sqrt((x - xc)**2 + (z - zc)**2) - R


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


    xp = xc[:,0,0]
    zp = zc[0,0,:]

    X, Z = np.meshgrid(xp, zp, indexing="ij")

    Sth = theta[:,0,:,j]
    Sth = np.ma.masked_where((Sth < 0) | np.isnan(Sth), Sth)

    fig, ax = plt.subplots(figsize=(7,3.5))

    levels = np.linspace(280, 320, 200)

    plot1 = ax.contourf(X, Z, Sth, levels=levels, cmap='RdBu_r', vmin=280, vmax=320)
    ax.contour(X, Z, Sth, levels=np.linspace(280, 320, 16), colors="k",linewidths=0.5)

    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_xlim([0, SizeX])
    ax.set_ylim([0, SizeZ])
    ax.set_aspect('equal', 'box')

    # ---- Dibujar círculo ----
    theta_circ = np.linspace(0, 2*np.pi, 300)
    x_circ = xc_circ + R_circ*np.cos(theta_circ)
    z_circ = zc_circ + R_circ*np.sin(theta_circ)


    cbar = plt.colorbar(plot1)
    cbar.ax.set_title('θ (K)')

    ax.fill(x_circ, z_circ, color='white', zorder=20)
    ax.plot(x_circ, z_circ, 'k', linewidth=1.0, zorder=21)


    plt.tight_layout()

    fig.savefig(fname+"_theta.png", dpi=300)
    plt.close(fig)  


    # =====================================================
    # Third plot: Zoom around circle with visible cells
    # =====================================================

    # ---- Define automatic zoom box around circle ----
    margin = 500.0

    x1_zoom = xc_circ - R_circ - margin
    x2_zoom = xc_circ + R_circ + margin
    z1_zoom = 0.0
    z2_zoom = 10000

    # ---- Compute cell edges from centers ----
    dx = xp[1] - xp[0]
    dz = zp[1] - zp[0]

    x_edges = np.concatenate(([xp[0] - dx/2], xp + dx/2))
    z_edges = np.concatenate(([zp[0] - dz/2], zp + dz/2))

    X_edges, Z_edges = np.meshgrid(x_edges, z_edges, indexing="ij")

    # ---- Extract indices inside zoom region ----
    ix = np.where((xp >= x1_zoom) & (xp <= x2_zoom))[0]
    iz = np.where((zp >= z1_zoom) & (zp <= z2_zoom))[0]

    ix_edges = np.arange(ix[0], ix[-1] + 2)
    iz_edges = np.arange(iz[0], iz[-1] + 2)

    Sth_zoom = Sth[ix[0]:ix[-1]+1, iz[0]:iz[-1]+1]

    X_edges_zoom = X_edges[ix_edges[0]:ix_edges[-1]+1,
                           iz_edges[0]:iz_edges[-1]+1]

    Z_edges_zoom = Z_edges[ix_edges[0]:ix_edges[-1]+1,
                           iz_edges[0]:iz_edges[-1]+1]

    Xp_zoom, Zp_zoom = np.meshgrid(
        xp[ix[0]:ix[-1]+1],
        zp[iz[0]:iz[-1]+1],
        indexing="ij"
    )

    # ---- Create figure ----
    fig3, ax3 = plt.subplots(figsize=(6,5))

    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_under('gray')

    plot3 = ax3.pcolormesh(
        X_edges_zoom,
        Z_edges_zoom,
        Sth_zoom,
        cmap=cmap,
        shading='flat',
        edgecolors='k',
        linewidth=0.2,
        vmin=280,
        vmax=320
    )

    # =====================================================
    # Mask centers using circle geometry (SDF)
    # =====================================================

    mask_circle = (Xp_zoom - xc_circ)**2 + (Zp_zoom - zc_circ)**2 <= R_circ**2

    # ---- Mask valid theta ----
    theta_mask = (Sth_zoom > 0) & (~np.isnan(Sth_zoom))

    mask_centers = mask_circle & theta_mask

    # ---- Plot cell centers inside circle ----
    ax3.scatter(
        Xp_zoom[mask_centers],
        Zp_zoom[mask_centers],
        s=6,
        marker='o',
        facecolors='none',
        edgecolors='black',
        linewidths=0.3,
        zorder=35
    )

    # ---- Gray cells (invalid θ or inside solid) ----
    gray_mask = ~theta_mask

    ax3.scatter(
        Xp_zoom[gray_mask],
        Zp_zoom[gray_mask],
        s=6,
        marker='o',
        color='black',
        linewidths=0.3,
        zorder=40
    )

    # =====================================================
    # Draw circle boundary
    # =====================================================

    theta_plot = np.linspace(0, 2*np.pi, 300)
    x_circ = xc_circ + R_circ*np.cos(theta_plot)
    z_circ = zc_circ + R_circ*np.sin(theta_plot)

    ax3.plot(x_circ, z_circ, 'k-', linewidth=1.0)

    # ---- Axis formatting ----
    ax3.set_xlim([x1_zoom, x2_zoom])
    ax3.set_ylim([z1_zoom, z2_zoom])
    ax3.set_aspect('equal', 'box')
    ax3.set_xlabel("x (m)")
    ax3.set_ylabel("z (m)")

    cbar3 = plt.colorbar(plot3)
    cbar3.ax.set_title('θ (K)')

    plt.tight_layout()

    zoom_name = fname + "_theta_cells"
    zoom_path = zoom_name + ".png"

    fig3.savefig(zoom_path, dpi=300)
    plt.close(fig3)    
    
    j=j+1
         
gif_path = os.path.join(folder_out, "animation.gif")
imageio.mimsave(gif_path, images, duration=8, loop=0)  # Adjust the duration as needed
