import math                    
import numpy as np             
import matplotlib.pyplot as plt 

folder_case="case"
folder_out="out"

fname_config="configure.input"
fname_eq="equilibrium.out"
fname_ini="initial.out"

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
Face_1 = 1 #-y
Face_2 = 4 #+x
Face_3 = 1 #+y
Face_4 = 4 #-x
Face_5 = 4 #-z
Face_6 = 4 #+z

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
v=np.zeros((xcells,ycells,zcells))
w=np.zeros((xcells,ycells,zcells))
rho=np.zeros((xcells,ycells,zcells))
p=np.zeros((xcells,ycells,zcells))  
phi=np.zeros((xcells,ycells,zcells)) 

ue=np.zeros((xcells,ycells,zcells))
ve=np.zeros((xcells,ycells,zcells))
we=np.zeros((xcells,ycells,zcells))
rhoe=np.zeros((xcells,ycells,zcells))
pe=np.zeros((xcells,ycells,zcells))  

    
x=np.arange(0+dx/2.0, SizeX, dx)
y=np.arange(0+dy/2.0, SizeY, dy)
z=np.arange(0+dz/2.0, SizeZ, dz)

xc, yc, zc= np.meshgrid(x,y,z,indexing='ij')

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

                u_i=0.0
                v_i=0.0
                w_i=0.0
                phi_i=0.0
                p_i=   p0*  (1.0-aux2*zc[l,m,n])**(gamma/(gamma-1.0))
                rho_i= rho0*(1.0-aux2*zc[l,m,n])**(1.0/(gamma-1.0))
                
                rhoe[l,m,n]=rho_i
                pe  [l,m,n]=p_i
                ue  [l,m,n]=u_i
                ve  [l,m,n]=v_i
                we  [l,m,n]=w_i

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
                p_i=   p0* (1.0-aux2*zc[l,m,n])**(gamma/(gamma-1.0)); 
                rho_i= p0/(R*tt)*(1.0-aux2*zc[l,m,n])**(1.0/(gamma-1.0));
                u_i=0.0;
                v_i=0.0;
                w_i=0.0;
                phi_i=0.0;
                
                rho[l,m,n]=rho_i
                p  [l,m,n]=p_i
                u  [l,m,n]=u_i
                v  [l,m,n]=v_i
                w  [l,m,n]=w_i
                phi[l,m,n]=phi_i


   
                
f = open(fname_config, "w")
f.write("/////SIMULATION_SETUP////// \n")
f.write("FinalTime    "+str(FinalTime)+"\n")
f.write("DumpTime    "+str(DumpTime)+"\n")
f.write("CFL    "+str(CFL)+"\n")
f.write("Order    "+str(Order)+"\n")
f.write(" \n")
f.write("////////MESH_SETUP/////////\n")
f.write("xcells    "+str(xcells)+"\n")
f.write("ycells    "+str(ycells)+"\n")
f.write("zcells    "+str(zcells)+"\n")
f.write("SizeX    "+str(SizeX)+"\n")
f.write("SizeY    "+str(SizeY)+"\n")
f.write("SizeZ    "+str(SizeZ)+"\n")
f.write(" \n")
f.write("///////BOUNDARY_COND///////\n")
f.write("Face_1    "+str(Face_1)+"\n")
f.write("Face_2    "+str(Face_2)+"\n")
f.write("Face_3    "+str(Face_3)+"\n")
f.write("Face_4    "+str(Face_4)+"\n")
f.write("Face_5    "+str(Face_5)+"\n")
f.write("Face_6    "+str(Face_6)+"\n")
f.write(" \n")
f.write("///////LINEAR_TRANSPORT///////(if_applicable)\n")
f.write("u_x    "+str(u_x)+"\n")
f.write("u_y    "+str(u_y)+"\n")
f.write("u_z    "+str(u_z)+"\n")
f.close()


f = open(fname_eq, "w")
f.write("VARIABLES = X, Y, Z, ue, ve, we, rhoe, pe, phi(n/u) \n")
f.write("CELLS = "+str(xcells)+", "+str(ycells)+", "+str(zcells)+","+"\n")
for l in range(0,xcells):   
    for m in range(0,ycells):
        for n in range(0,zcells):
            f.write(str(xc[l,m,n])+" "+str(yc[l,m,n])+" "+str(zc[l,m,n])+" "+str(ue[l,m,n])+" "+str(ve[l,m,n])+" "+str(we[l,m,n])+" "+str(rhoe[l,m,n])+" "+str(pe[l,m,n])+" 0.0"+"\n")

f.close()         




f = open(fname_ini, "w")
f.write("VARIABLES = X, Y, Z, u, v, w, rho, p, phi \n")
f.write("CELLS = "+str(xcells)+", "+str(ycells)+", "+str(zcells)+","+"\n")
for l in range(0,xcells):   
    for m in range(0,ycells):
        for n in range(0,zcells):
            f.write(str(xc[l,m,n])+" "+str(yc[l,m,n])+" "+str(zc[l,m,n])+" "+str(u[l,m,n])+" "+str(v[l,m,n])+" "+str(w[l,m,n])+" "+str(rho[l,m,n])+" "+str(p[l,m,n])+" "+str(phi[l,m,n])+"\n")

f.close()         
