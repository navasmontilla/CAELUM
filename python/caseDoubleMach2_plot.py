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
folder_case="caseDoubleMach/" 


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




#Mesh setup
xcells = 300
ycells = 200
zcells = 1
SizeX = 3.0
SizeY = 2.0
SizeZ = SizeY/ycells*zcells
xc, yc, zc, u, v, w, rho, p, phi, ue, ve, we, rhoe, pe = initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ)



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



    # =====================================================
    # Third plot: Zoomed square cells with visible edges
    # =====================================================

    # ---- Define zoom box ----
    x1, x2 = 2.4, 2.8
    y1, y2 = 1.0, 1.4

    # ---- Compute cell edges from centers (uniform grid assumed) ----
    dx = xp[1] - xp[0]
    dy = yp[1] - yp[0]

    x_edges = np.concatenate(([xp[0] - dx/2], xp + dx/2))
    y_edges = np.concatenate(([yp[0] - dy/2], yp + dy/2))

    X_edges, Y_edges = np.meshgrid(x_edges, y_edges)

    # ---- Extract indices inside zoom region ----
    ix = np.where((xp >= x1) & (xp <= x2))[0]
    iy = np.where((yp >= y1) & (yp <= y2))[0]

    # edges need +1 in upper bound
    ix_edges = np.arange(ix[0], ix[-1] + 2)
    iy_edges = np.arange(iy[0], iy[-1] + 2)

    Srho_zoom = Srho[iy[0]:iy[-1]+1, ix[0]:ix[-1]+1]

    X_edges_zoom = X_edges[iy_edges[0]:iy_edges[-1]+1,
                            ix_edges[0]:ix_edges[-1]+1]
    Y_edges_zoom = Y_edges[iy_edges[0]:iy_edges[-1]+1,
                            ix_edges[0]:ix_edges[-1]+1]

    # ---- Create center grid for markers ----
    Xp_zoom, Yp_zoom = np.meshgrid(
        xp[ix[0]:ix[-1]+1],
        yp[iy[0]:iy[-1]+1]
    )

    # ---- Create figure ----
    fig3, ax3 = plt.subplots(figsize=(5,4))

    # Colormap with gray for rho < 0
    cmap = plt.get_cmap('RdBu').copy()
    cmap.set_under('gray')

    plot3 = ax3.pcolormesh(
        X_edges_zoom,
        Y_edges_zoom,
        Srho_zoom,
        cmap=cmap,
        shading='flat',
        edgecolors='k',
        linewidth=0.2,
        vmin=0.0,
        vmax=14
    )

    # ---- 30 degree line ----
    m = np.tan(np.deg2rad(30))
    x_line_zoom = np.linspace(x1, x2, 500)
    y_line_zoom = m * (x_line_zoom - 0.5)

    ax3.plot(x_line_zoom, y_line_zoom, color='black', linewidth=0.6)

    # ---- Mask for markers ----
    line_mask = Yp_zoom < m*(Xp_zoom - 0.5)
    rho_mask = Srho_zoom >= 0
    mask_centers = line_mask & rho_mask

    # ---- Plot cell center markers ----
    ax3.scatter(
        Xp_zoom[mask_centers],
        Yp_zoom[mask_centers],
        s=6,
        marker='o',
        facecolors='none',
        edgecolors='black',
        linewidths=0.3,
        zorder=20
    )
    
    # ---- Mask for gray cells (rho < 0) ----
    gray_mask = Srho_zoom < 0

    # ---- Plot crosses on gray cells ----
    ax3.scatter(
        Xp_zoom[gray_mask],
        Yp_zoom[gray_mask],
        s=6,
        marker='o',
        color='black',
        linewidths=0.3,
        zorder=25
)

    # ---- Axis formatting ----
    ax3.set_xlim([x1, x2])
    ax3.set_ylim([y1, y2])
    ax3.set_aspect('equal', 'box')
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")

    cbar3 = plt.colorbar(plot3)
    cbar3.ax.set_title('ρ')    
    

    plt.tight_layout()

    zoom_name = fname + "_doublemach_cells"
    zoom_path = zoom_name + ".png"

    fig3.savefig(zoom_path, dpi=300)
    plt.close(fig3)
    
    
    j=j+1
         

gif_path = os.path.join(folder_out, "animation.gif")
imageio.mimsave(gif_path, images, duration=0.01, loop=0)  # Adjust the duration as needed