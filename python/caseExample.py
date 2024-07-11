import os
import math                    
import numpy as np             
import matplotlib.pyplot as plt 
import re
from glob import glob
from utils import modify_header_file,write_config,write_initial,write_equilibrium,backup_file,restore_file,compile_program,run_program,initialize_variables,read_data_euler

script_dir = os.path.dirname(os.path.abspath(__file__))

folder_case="caseExample/" 

#Do not modify the folders and paths below
folder_out="out/"
folder_lib="lib"

fname_config="configure.input"
fname_eq="equilibrium.out"
fname_ini="initial.out"

folder_out = os.path.join(script_dir, "../"+folder_case+"/"+folder_out)
folder_case = os.path.join(script_dir, "../"+folder_case)
folder_lib = os.path.join(script_dir, "../"+folder_lib)
folder_exe = os.path.join(script_dir, "../")

backup_file(folder_lib+'/definitions.h')

#Configure the header file for compilation
modify_header_file(folder_lib+'/definitions.h', 'NTRHEADS', 2)
modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 2)
modify_header_file(folder_lib+'/definitions.h', 'ST', 3)
modify_header_file(folder_lib+'/definitions.h', 'SOLVER', 0)
modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)

#Simulation setup
FinalTime = 800.0
DumpTime = 100.0
CFL = 0.45
Order = 7

#Mesh setup
xcells = 100
ycells = 1
zcells = 50
SizeX = 20000.0
SizeY = 1000.0
SizeZ = 10000.0

#Boundary conditions
Face_1 = 4 #-y
Face_2 = 4 #+x
Face_3 = 4 #+y
Face_4 = 4 #-x
Face_5 = 4 #-z
Face_6 = 4 #+z

#Linear transport, only if applicable
u_x = 1.0
u_y = 1.0
u_z = 1.0


xc, yc, zc, u, v, w, rho, p, phi, ue, ve, we, rhoe, pe = initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ)


#INITIAL CONDITION

tt0=300
p0=1.e5
R=287.058
gamma=1.4
g=9.8
rho0=p0/(R*tt0)
aux2=(gamma-1.0)/gamma*g/(R*tt0)
            
for l in range(0,xcells): 
        for m in range(0,ycells): 
            for n in range(0,zcells):
                
                rhoe[l,m,n]=rho0*(1.0-aux2*zc[l,m,n])**(1.0/(gamma-1.0))
                pe  [l,m,n]=p0*  (1.0-aux2*zc[l,m,n])**(gamma/(gamma-1.0))
                ue  [l,m,n]=0.0
                ve  [l,m,n]=0.0
                we  [l,m,n]=0.0

                xp=10000;
                zp=2000;
                d1=np.sqrt((xc[l,m,n]-xp)*(xc[l,m,n]-xp)+(zc[l,m,n]-zp)*(zc[l,m,n]-zp));

                xp=10000;
                zp=8000;
                d2=np.sqrt((xc[l,m,n]-xp)*(xc[l,m,n]-xp)+(zc[l,m,n]-zp)*(zc[l,m,n]-zp));

                rc=1000;
                aux1=20.0*(max(rc-d1/2.0,0.0)+min(d2/2.0-rc,0.0))/1000;

                tt=tt0+aux1;
                aux2=(gamma-1.0)/gamma*g/(R*tt0);
                
                rho[l,m,n]=p0/(R*tt)*(1.0-aux2*zc[l,m,n])**(1.0/(gamma-1.0))
                p  [l,m,n]=p0* (1.0-aux2*zc[l,m,n])**(gamma/(gamma-1.0))
                u  [l,m,n]=0.0
                v  [l,m,n]=0.0
                w  [l,m,n]=0.0
                phi[l,m,n]=0.0

  
# WRITING CONFIGURATION AND INITIAL DATA

write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, Order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z)
write_equilibrium(folder_case, fname_eq, xcells, ycells, zcells, xc, yc, zc, ue, ve, we, rhoe, pe)       
write_initial(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u, v, w, rho, p, phi)

# COMPILING AND RUNNING

compile_program()
print("Program is running...")
run_program(folder_exe+"./exehow3d "+folder_case)

# READING OUTPUT DATA AND PLOTTING

files = glob(folder_out+"/*.out")
lf=len(files)
print(files)

gamma=1.4
j=0
print("Printing figures in folder"+folder_out)
for fname in files:
    
    u, v, w, rho, p, phi, theta, E = read_data_euler(fname, xcells, ycells, zcells, lf, gamma, j)
                  
    filename = fname+"_theta2D"
    
    xp = xc[:,0,0]     #puntos en x
    zp = zc[0,0,:]      #puntos en y
    X, Y = np.meshgrid(xp, zp)    #matriz de puntos
    Srho=np.transpose(rho[:,0,:,j])
    Spres=np.transpose(p[:,0,:,j])
    Svel=np.transpose(np.sqrt(u[:,0,:,j]**2+w[:,0,:,j]**2))
    Senr=np.transpose(E[:,0,:,j])
    Sphi=np.transpose(phi[:,0,:,j])
    Su=np.transpose(u[:,0,:,j])
    Sth=np.transpose(theta[:,0,:,j])
    
    fig, ax = plt.subplots(figsize=(10, 5))      #genera el objeto "figura"
    levels = np.linspace(280, 320, 16)
    #print(levels)
    plot1=ax.contour(X, Y, Sth, levels=levels,colors="k",linewidths=0.5)  
    plot1=ax.contourf(X, Y, Sth, 200, cmap='RdBu')   
    ax.set_title('Potential temperature')
    ax.set_xlabel("x") 
    ax.set_ylabel("z") 
    ax.set_aspect('equal', 'box')
    plot1.set_clim( 280, 320 )
    # Create colorbar
    cbar = plt.colorbar(plot1)
    #cbar_ticks = levels
    #cbar.set_ticks(cbar_ticks)
    cbar.ax.set_title('θ(K)')
    fig.text(0.15, 0.72, "θ_max="+str(round(np.max(Sth),3)), fontsize=9.5)
    fig.text(0.15, 0.68, "θ_min="+str(round(np.min(Sth),3)), fontsize=9.5)
    
    fig.savefig(filename+".png",dpi=500)
    
    j=j+1
         

restore_file(folder_lib+'/definitions.h')