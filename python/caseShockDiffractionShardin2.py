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
from matplotlib.patches import Polygon
import re
from glob import glob
from utils import modify_header_file,write_config,write_initial,write_equilibrium,backup_file,restore_file,compile_program,run_program,initialize_variables,read_data_euler,write_solid_cells
import imageio
import pyvista as pv

def udTriangle2D(p, a, b, c):
    """
    Signed distance function to triangle in 2D.
    Negative inside, positive outside.
    p = np.array([x, y])
    a, b, c = vertices of triangle (np.array([x,y]))
    """
    # Distancia mínima a los lados
    def seg_sdf(pv, v0, v1):
        v = v1 - v0
        w = pv - v0
        t = np.clip(np.dot(w, v)/np.dot(v,v), 0.0, 1.0)
        proj = v0 + t*v
        return np.linalg.norm(pv - proj)

    d = min(seg_sdf(p, a, b), seg_sdf(p, b, c), seg_sdf(p, c, a))

    # Coordenadas baricéntricas
    v0 = b - a
    v1 = c - a
    v2 = p - a

    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)

    invDenom = 1.0 / (dot00*dot11 - dot01*dot01)
    u = (dot11*dot02 - dot01*dot12) * invDenom
    v = (dot00*dot12 - dot01*dot02) * invDenom

    if (u >= 0) and (v >= 0) and (u + v <= 1):
        return -d  # negativo dentro
    else:
        return d   # positivo fuera


# ### Setting up the paths
# 
# First, the name of the folder for this test case must be specified:

#This test case will run in the folder "caseRP3D/". 
#Don't forget the bar (/). 
#This directory should have been created prior to the execution to this script, and should also contain an /out folder inside
folder_case="caseShockDiffractionSchardin2/" 


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
modify_header_file(folder_lib+'/definitions.h', 'SOLVER', 1)           #Riemann solver used
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
FinalTime = 0.000151
DumpTime = 0.00005
CFL = 0.4 #0.4 in 2D
Order = 7

#Mesh setup
xcells = 2000
ycells = 1600
zcells = 1
SizeX = 0.15
SizeY = 0.12
SizeZ = SizeY/ycells*zcells

#Boundary conditions
Face_1 = 3 #-y
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

# -------------------------------------------------
# Physical parameters (Schardin case)
# -------------------------------------------------
gamma = 1.4
Ms = 1.34
pR = 0.05e6
TR = 285.0
Rgas = 287.0

rhoR = pR/(Rgas*TR)
uR = 0.0
vR = 0.0

aR = np.sqrt(gamma*pR/rhoR)
shock_speed = Ms*aR

rhoL = rhoR*((gamma+1)*Ms**2)/((gamma-1)*Ms**2 + 2)
pL = pR*(1 + 2*gamma/(gamma+1)*(Ms**2 - 1))

uL = shock_speed*(1 - rhoR/rhoL)
vL = 0.0

# -------------------------------------------------
# Shock initial position
# -------------------------------------------------
xs = 0.050   # 2 cm from inlet (adjust if desired)

# -------------------------------------------------
# Triangle geometry (30° half-angle)
# -------------------------------------------------
theta = np.deg2rad(60.0/2.0)

xa = 0.050          # cusp location (m)
ya = SizeY/2.0      # centered vertically
length_x = 0.010/np.tan(theta)    #  cm triangle length

slope = np.tan(theta)

x_base = xa + length_x   # ahora la base está a la derecha

# Triangle vertices
x0, y0 = xa, ya                  # vertex (cusp)
x1, y1 = x_base, ya + slope*length_x   # top base
x2, y2 = x_base, ya - slope*length_x   # bottom base

# Triangle vertices
v0 = np.array([x0, y0])
v1 = np.array([x1, y1])
v2 = np.array([x2, y2])

for l in range(xcells):
    for m in range(ycells):
        for n in range(zcells):

            x = xc[l,m,n]
            y = yc[l,m,n]
            p_point = np.array([x, y])

            # -------------------------
            # Planar moving shock
            # -------------------------
            if x < xs:
                p[l,m,n]   = pL
                rho[l,m,n] = rhoL
                u[l,m,n]   = uL
                v[l,m,n]   = vL
                phi[l,m,n] = 1.0
            else:
                p[l,m,n]   = pR
                rho[l,m,n] = rhoR
                u[l,m,n]   = uR
                v[l,m,n]   = vR
                phi[l,m,n] = 0.0

            # -------------------------
            # Signed distance function
            # -------------------------
            sld[l,m,n] = udTriangle2D(p_point, v0, v1, v2)

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
                  
    filename = fname+"_schardin"
    
    xp = xc[:,0,0]     
    yp = yc[0,:,0]      
    X, Y = np.meshgrid(xp, yp)    #matriz de puntos
    Srho=np.transpose(rho[:,:,0,j])
    Spres=np.transpose(p[:,:,0,j])
    Svel=np.transpose(np.sqrt(u[:,:,0,j]**2+w[:,:,0,j]**2))
    Senr=np.transpose(E[:,:,0,j])
    Sphi=np.transpose(phi[:,:,0,j])
    Su=np.transpose(u[:,:,0,j])
    
    lev=np.linspace(0.1,1.1,512)
    fig, ax = plt.subplots(figsize=(7, 4))   
    #levels = np.linspace(0, 4, 10)
    #print(levels)
    #plot1=ax.contour(X, Y, Srho, levels=lev, colors="k",linewidths=0.2)  
    plot1=ax.contourf(X, Y, Srho, levels=lev, cmap='RdBu_r')   
    ax.set_xlabel("x") 
    ax.set_ylabel("y") 
    #ax.set_xlim([0,3]) 
    #ax.set_ylim([0,2]) 
    ax.set_aspect('equal', 'box')
    #plot1.set_clim( 2, 23 )
    # Create colorbar
    cbar = plt.colorbar(plot1)
    cbar.ax.set_title('ρ')
    plt.tight_layout()
    
    ddx=.08
    mm=57.9/80.0
    #ax.plot([x0,x0+ddx],[y0,y0+mm*ddx],'k',linewidth=0.2)
    
    # ---- Dibujar triángulo sólido ----
    triangle_coords = np.array([
        [x0, y0],
        [x1, y1],
        [x2, y2]
    ])

    triangle = Polygon(triangle_coords,
                       closed=True,
                       facecolor='white',
                       edgecolor='black',
                       linewidth=0.8,
                       zorder=20)

    ax.add_patch(triangle)
    
    image_path = filename + ".png"
    fig.savefig(image_path,dpi=250)
    images.append(imageio.imread(image_path)) 
    
    plt.close(fig)
   

    # =====================================================
    # Second plot: Density gradient magnitude |∇rho|
    # =====================================================

    # Compute spatial steps
    dx = xp[1] - xp[0]
    dy = yp[1] - yp[0]

    # Compute gradients (note: Srho is already transposed as [y,x])
    drho_dy, drho_dx = np.gradient(Srho, dy, dx)

    grad_rho = np.sqrt(drho_dx**2 + drho_dy**2)

    fig2, ax2 = plt.subplots(figsize=(7,6))

    levels_grad = np.linspace(0.0, 250, 512)

    cmap = plt.get_cmap('gray_r').copy()
    cmap.set_over('black')   # values > vmax will be black
    cmap.set_under('white')  # values < vmin will be white (optional)
    
    plot2 = ax2.imshow(
        grad_rho,
        extent=[xp.min(), xp.max(), yp.min(), yp.max()],
        origin='lower',
        cmap=cmap,
        vmin=0,
        vmax=250,
        aspect='equal'
    )

    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_aspect('equal', 'box')

    #cbar2 = plt.colorbar(plot2)
    #cbar2.ax.set_title('|∇ρ|')

    # Draw triangle (solid body)
    triangle_coords = np.array([
        [x0, y0],
        [x1, y1],
        [x2, y2]
    ])

    triangle = Polygon(
        triangle_coords,
        closed=True,
        facecolor='white',
        edgecolor='black',
        linewidth=0.8,
        zorder=20
    )

    ax2.add_patch(triangle)

    plt.tight_layout()

    grad_name = fname + "_grad_rho"
    grad_path = grad_name + ".png"

    fig2.savefig(grad_path, dpi=300)
    plt.close(fig2)   
    
    
    
    # =====================================================
    # Third plot: Zoom around triangle with visible cells
    # =====================================================

    # ---- Define automatic zoom box around triangle ----
    margin = 0.005

    x_min_tri = min(x0, x1, x2)
    x_max_tri = max(x0, x1, x2)
    y_min_tri = min(y0, y1, y2)
    y_max_tri = max(y0, y1, y2)

    x1_zoom = x_min_tri - margin
    x2_zoom = x_max_tri + margin
    y1_zoom = y_min_tri - margin
    y2_zoom = y_max_tri + margin

    # ---- Compute cell edges from centers ----
    dx = xp[1] - xp[0]
    dy = yp[1] - yp[0]

    x_edges = np.concatenate(([xp[0] - dx/2], xp + dx/2))
    y_edges = np.concatenate(([yp[0] - dy/2], yp + dy/2))

    X_edges, Y_edges = np.meshgrid(x_edges, y_edges)

    # ---- Extract indices inside zoom region ----
    ix = np.where((xp >= x1_zoom) & (xp <= x2_zoom))[0]
    iy = np.where((yp >= y1_zoom) & (yp <= y2_zoom))[0]

    ix_edges = np.arange(ix[0], ix[-1] + 2)
    iy_edges = np.arange(iy[0], iy[-1] + 2)

    Srho_zoom = Srho[iy[0]:iy[-1]+1, ix[0]:ix[-1]+1]

    X_edges_zoom = X_edges[iy_edges[0]:iy_edges[-1]+1,
                            ix_edges[0]:ix_edges[-1]+1]
    Y_edges_zoom = Y_edges[iy_edges[0]:iy_edges[-1]+1,
                            ix_edges[0]:ix_edges[-1]+1]

    Xp_zoom, Yp_zoom = np.meshgrid(
        xp[ix[0]:ix[-1]+1],
        yp[iy[0]:iy[-1]+1]
    )

    # ---- Create figure ----
    fig3, ax3 = plt.subplots(figsize=(5,4))

    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_under('gray')

    plot3 = ax3.pcolormesh(
        X_edges_zoom,
        Y_edges_zoom,
        Srho_zoom,
        cmap=cmap,
        shading='flat',
        edgecolors='k',
        linewidth=0.2,
        vmin=0.1,
        vmax=1.1
    )

    # =====================================================
    # Draw triangle (filled white)
    # =====================================================

    triangle_coords = np.array([
        [x0, y0],
        [x1, y1],
        [x2, y2]
    ])

    triangle = Polygon(
        triangle_coords,
        closed=True,
        facecolor='none',
        edgecolor='black',
        linewidth=1.0,
        zorder=30
    )

    #ax3.add_patch(triangle)

    # =====================================================
    # Mask centers using triangle geometry
    # =====================================================

    def point_in_triangle(px, py, a, b, c):
        v0 = c - a
        v1 = b - a
        v2 = np.array([px, py]) - a

        dot00 = np.dot(v0, v0)
        dot01 = np.dot(v0, v1)
        dot02 = np.dot(v0, v2)
        dot11 = np.dot(v1, v1)
        dot12 = np.dot(v1, v2)

        invDen = 1.0 / (dot00*dot11 - dot01*dot01)
        u = (dot11*dot02 - dot01*dot12) * invDen
        v = (dot00*dot12 - dot01*dot02) * invDen

        return (u >= 0) & (v >= 0) & (u + v <= 1)

    v0 = np.array([x0, y0])
    v1 = np.array([x1, y1])
    v2 = np.array([x2, y2])

    mask_triangle = np.zeros_like(Srho_zoom, dtype=bool)

    for i in range(Xp_zoom.shape[0]):
        for j2 in range(Xp_zoom.shape[1]):
            mask_triangle[i, j2] = point_in_triangle(
                Xp_zoom[i, j2],
                Yp_zoom[i, j2],
                v0, v1, v2
            )

    rho_mask = Srho_zoom >= 0
    mask_centers = mask_triangle & rho_mask

    # ---- Plot cell centers inside triangle ----
    ax3.scatter(
        Xp_zoom[mask_centers],
        Yp_zoom[mask_centers],
        s=6,
        marker='o',
        facecolors='none',
        edgecolors='black',
        linewidths=0.3,
        zorder=35
    )

    # ---- Gray cells (rho < 0) ----
    gray_mask = Srho_zoom < 0

    ax3.scatter(
        Xp_zoom[gray_mask],
        Yp_zoom[gray_mask],
        s=6,
        marker='o',
        color='black',
        linewidths=0.3,
        zorder=40
    )

    # ---- Axis formatting ----
    ax3.set_xlim([x1_zoom, x2_zoom])
    ax3.set_ylim([y1_zoom, y2_zoom])
    ax3.set_aspect('equal', 'box')
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")

    cbar3 = plt.colorbar(plot3)
    cbar3.ax.set_title('ρ')

    plt.tight_layout()

    zoom_name = fname + "_schardin_cells"
    zoom_path = zoom_name + ".png"

    fig3.savefig(zoom_path, dpi=300)
    plt.close(fig3)
    
    
    
    
    
    
    
    j=j+1
         

gif_path = os.path.join(folder_out, "animation.gif")
imageio.mimsave(gif_path, images, duration=8, loop=0)  # Adjust the duration as needed