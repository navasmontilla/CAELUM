import math                    
import numpy as np             
import matplotlib.pyplot as plt 

folder_case="case"
folder_out="out"

fname_config="configure.input"
fname_ini="initial.out"

#Simulation setup
FinalTime = 0.4
DumpTime = 0.9
CFL = 0.3
Order = 5

#Mesh setup
xcells = 40
ycells = 40
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


dx=SizeX/xcells
dy=SizeY/ycells
dz=SizeZ/zcells

xc=np.zeros((xcells,ycells,zcells))
yc=np.zeros((xcells,ycells,zcells))
zc=np.zeros((xcells,ycells,zcells))

u=np.zeros((xcells,ycells,zcells)) 
    
x=np.arange(0+dx/2.0, SizeX, dx)
y=np.arange(0+dy/2.0, SizeY, dy)
z=np.arange(0+dz/2.0, SizeZ, dz)

xc, yc, zc= np.meshgrid(x,y,z,indexing='ij')

for l in range(0,xcells): 
        for m in range(0,ycells):
            r=np.sqrt((cell[k].xc-xc)*(cell[k].xc-xc)+(cell[k].yc-yc)*(cell[k].yc-yc));
            if (r<0.5):
                rho[l,m,:]=77.0/558.0
                p  [l,m,:]=9.0/310.0
                u  [l,m,:]=4.0/np.sqrt(11.0)
                v  [l,m,:]=4.0/np.sqrt(11.0)
                phi[l,m,:]=1.0