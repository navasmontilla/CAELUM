# utils.py

import os
import re
import subprocess
import shutil
import math
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
from glob import glob

def modify_header_file(file_path, macro_name, new_value):
    """
    Modifies the definition of a macro in a header file.
    
    Args:
        file_path (str): Path to the header file.
        macro_name (str): Name of the macro to be modified.
        new_value (str/int/float): New value for the macro.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    modified_lines = []
    for line in lines:
        if re.match(rf"#define\s+{macro_name}\s+", line):
            line = re.sub(rf"#define\s+{macro_name}\s+\d+", f"#define {macro_name} {new_value}", line)
        modified_lines.append(line)
    
    with open(file_path, 'w') as file:
        file.writelines(modified_lines)
    print(f"{macro_name} changed to {new_value} in {file_path}")


def backup_file(file_path):
    """
    Creates a backup of the specified file.
    
    Args:
        file_path (str): Path to the file to be backed up.
    """
    backup_path = file_path + '.bak'
    shutil.copy(file_path, backup_path)
    print(f"Backup created at {backup_path}")

def restore_file(file_path):
    """
    Restores the original file from its backup and removes the backup.
    
    Args:
        file_path (str): Path to the file to be restored.
    """
    backup_path = file_path + '.bak'
    if os.path.exists(backup_path):
        shutil.copy(backup_path, file_path)
        os.remove(backup_path)
        print(f"File restored from {backup_path} and backup removed")
    else:
        print(f"No backup found at {backup_path}")
        
def backup_fileL2(file_path):
    """
    Creates a backup of the specified file.
    
    Args:
        file_path (str): Path to the file to be backed up.
    """
    backup_path = file_path + 'L2.bak'
    shutil.copy(file_path, backup_path)
    print(f"Backup created at {backup_path}")

def restore_fileL2(file_path):
    """
    Restores the original file from its backup and removes the backup.
    
    Args:
        file_path (str): Path to the file to be restored.
    """
    backup_path = file_path + 'L2.bak'
    if os.path.exists(backup_path):
        shutil.copy(backup_path, file_path)
        os.remove(backup_path)
        print(f"File restored from {backup_path} and backup removed")
    else:
        print(f"No backup found at {backup_path}")
        
        
def compile_program(make_command='make'):
    """
    Compiles the program using the specified command after running `make clean`,
    and streams the output to the console in real-time.
    
    Args:
        make_command (str): Command to compile the program. Default is 'make'.
    """
    # Run the clean command
    clean_command = 'make clean'
    print("Running clean command...")
    clean_result = subprocess.run(clean_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(clean_result.stdout)

    if clean_result.returncode != 0:
        print(f"Clean error: {clean_result.stdout}")
        return

    # Run the compile command
    print("Running compile command...")
    compile_process = subprocess.Popen(make_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    
    # Stream the output to the console
    for line in iter(compile_process.stdout.readline, ''):
        print(line, end='')

    compile_process.stdout.close()
    compile_process.wait()

    if compile_process.returncode != 0:
        print(f"Compilation error with exit code {compile_process.returncode}")
    else:
        print("Compilation successful")
        
def compile_program_jupyter(make_command='make', makefile_directory='../'): 
    """
    Compiles the program using the specified command after running `make clean`,
    and streams the output to the console in real-time.
    
    Args:
        make_command (str): Command to compile the program. Default is 'make'.
        makefile_directory (str): Directory where the Makefile is located. Default is the current directory.
    """
    # Run the clean command (if required)
    clean_command = 'make clean'
    print("Running clean command...")
    clean_result = subprocess.run(clean_command, shell=True, cwd=makefile_directory, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(clean_result.stdout)

    if clean_result.returncode != 0:
        print(f"Clean error: {clean_result.stdout}")
        return

    # Run the compile command
    print("Running compile command...")
    compile_process = subprocess.Popen(make_command, shell=True, cwd=makefile_directory, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    
    # Stream the output to the console
    for line in iter(compile_process.stdout.readline, ''):
        print(line, end='')

    compile_process.stdout.close()
    compile_process.wait()

    if compile_process.returncode != 0:
        print(f"Compilation error with exit code {compile_process.returncode}")
    else:
        print("Compilation successful")        

def run_program(executable_path):
    """
    Runs the compiled program and captures the output to display on the screen.
    
    Args:
        executable_path (str): Full path to the program's executable.
    """
    result = subprocess.run(executable_path, shell=True, capture_output=True, text=True)
    
    # Print the captured output and error (if any)
    print("Program output:")
    print(result.stdout)
    
        
def initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ):

    dx=SizeX/xcells
    dy=SizeY/ycells
    dz=SizeZ/zcells

    xc = np.zeros((xcells, ycells, zcells))
    yc = np.zeros((xcells, ycells, zcells))
    zc = np.zeros((xcells, ycells, zcells))

    u = np.zeros((xcells, ycells, zcells))
    v = np.zeros((xcells, ycells, zcells))
    w = np.zeros((xcells, ycells, zcells))
    rho = np.zeros((xcells, ycells, zcells))
    p = np.zeros((xcells, ycells, zcells))
    phi = np.zeros((xcells, ycells, zcells))

    ue = np.zeros((xcells, ycells, zcells))
    ve = np.zeros((xcells, ycells, zcells))
    we = np.zeros((xcells, ycells, zcells))
    rhoe = np.zeros((xcells, ycells, zcells))
    pe = np.zeros((xcells, ycells, zcells))

    x = np.arange(0 + dx / 2.0, SizeX, dx)
    y = np.arange(0 + dy / 2.0, SizeY, dy)
    z = np.arange(0 + dz / 2.0, SizeZ, dz)

    xc, yc, zc = np.meshgrid(x, y, z, indexing='ij')

    return xc, yc, zc, u, v, w, rho, p, phi, ue, ve, we, rhoe, pe


def write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, Order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z):
    with open(folder_case + fname_config, "w") as f:
        f.write("/////SIMULATION_SETUP////// \n")
        f.write(f"FinalTime    {FinalTime}\n")
        f.write(f"DumpTime    {DumpTime}\n")
        f.write(f"CFL    {CFL}\n")
        f.write(f"Order    {Order}\n")
        f.write(" \n")
        f.write("////////MESH_SETUP/////////\n")
        f.write(f"xcells    {xcells}\n")
        f.write(f"ycells    {ycells}\n")
        f.write(f"zcells    {zcells}\n")
        f.write(f"SizeX    {SizeX}\n")
        f.write(f"SizeY    {SizeY}\n")
        f.write(f"SizeZ    {SizeZ}\n")
        f.write(" \n")
        f.write("///////BOUNDARY_COND///////\n")
        f.write(f"Face_1    {Face_1}\n")
        f.write(f"Face_2    {Face_2}\n")
        f.write(f"Face_3    {Face_3}\n")
        f.write(f"Face_4    {Face_4}\n")
        f.write(f"Face_5    {Face_5}\n")
        f.write(f"Face_6    {Face_6}\n")
        f.write(" \n")
        f.write("///////LINEAR_TRANSPORT///////(if_applicable)\n")
        f.write(f"u_x    {u_x}\n")
        f.write(f"u_y    {u_y}\n")
        f.write(f"u_z    {u_z}\n")


def write_equilibrium(folder_case, fname_eq, xcells, ycells, zcells, xc, yc, zc, ue, ve, we, rhoe, pe):
    with open(folder_case + fname_eq, "w") as f:
        f.write("VARIABLES = X, Y, Z, ue, ve, we, rhoe, pe, phi(n/u) \n")
        f.write(f"CELLS = {xcells}, {ycells}, {zcells},\n")
        for l in range(xcells):
            for m in range(ycells):
                for n in range(zcells):
                    f.write(f"{xc[l,m,n]} {yc[l,m,n]} {zc[l,m,n]} {ue[l,m,n]} {ve[l,m,n]} {we[l,m,n]} {rhoe[l,m,n]} {pe[l,m,n]} 0.0\n")


def write_initial(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u, v, w, rho, p, phi):
    with open(folder_case + fname_ini, "w") as f:
        f.write("VARIABLES = X, Y, Z, u, v, w, rho, p, phi \n")
        f.write(f"CELLS = {xcells}, {ycells}, {zcells},\n")
        for l in range(xcells):
            for m in range(ycells):
                for n in range(zcells):
                    f.write(f"{xc[l,m,n]} {yc[l,m,n]} {zc[l,m,n]} {u[l,m,n]} {v[l,m,n]} {w[l,m,n]} {rho[l,m,n]} {p[l,m,n]} {phi[l,m,n]}\n")
                    
def write_initial_scalar(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u):
    with open(folder_case + fname_ini, "w") as f:
        f.write("VARIABLES = X, Y, Z, u \n")
        f.write(f"CELLS = {xcells}, {ycells}, {zcells},\n")
        for l in range(xcells):
            for m in range(ycells):
                for n in range(zcells):
                    f.write(f"{xc[l,m,n]} {yc[l,m,n]} {zc[l,m,n]} {u[l,m,n]} \n")

def write_solid_cells(folder_case, fname_solids, xcells, ycells, zcells, xc, yc, zc, u):
    with open(folder_case + fname_solids, "w") as f:
        f.write("VARIABLES = X, Y, Z, u \n")
        f.write(f"CELLS = {xcells}, {ycells}, {zcells},\n")
        for l in range(xcells):
            for m in range(ycells):
                for n in range(zcells):
                    f.write(f"{xc[l,m,n]} {yc[l,m,n]} {zc[l,m,n]} {u[l,m,n]} \n")                    

def read_data_euler(fname, xcells, ycells, zcells, lf, gamma, j):
    """
    Process data from a file and return processed arrays.
    
    Args:
                
    Returns:
        tuple: Arrays (u, v, w, rho, p, phi, theta) processed from the data file.
    """
    u = np.zeros((xcells, ycells, zcells, lf))
    v = np.zeros((xcells, ycells, zcells, lf))
    w = np.zeros((xcells, ycells, zcells, lf))
    rho = np.zeros((xcells, ycells, zcells, lf))
    p = np.zeros((xcells, ycells, zcells, lf))
    phi = np.zeros((xcells, ycells, zcells, lf))
    theta = np.zeros((xcells, ycells, zcells, lf))
    E = np.zeros((xcells, ycells, zcells, lf))
    rhoE = np.zeros((xcells, ycells, zcells))
    pE = np.zeros((xcells, ycells, zcells))

    file1 = open(fname, 'r')
    print(fname + " file read")
    data = np.loadtxt(fname, skiprows=2)
    k = 0
    for l in range(0, xcells):
        for m in range(0, ycells):
            for n in range(0, zcells):
                u[l, m, n, j] = data[k, 3]
                v[l, m, n, j] = data[k, 4]
                w[l, m, n, j] = data[k, 5]
                rho[l, m, n, j] = data[k, 6]
                p[l, m, n, j] = data[k, 7]
                phi[l, m, n, j] = data[k, 8]
                theta[l, m, n, j] = data[k, 9]
                E[l,m,n,j]=p[l,m,n,j]/(gamma-1.0)+0.5*rho[l,m,n,j]*(u[l,m,n,j]*u[l,m,n,j] + v[l,m,n,j]*v[l,m,n,j] + w[l,m,n,j]*w[l,m,n,j]);
                k += 1

    return u, v, w, rho, p, phi, theta, E
    
def read_data_scalar(u,fname, xcells, ycells, zcells, lf, j):
    """
    Process data from a file and return processed arrays.
    
    Args:
                
    Returns:
        tuple: Arrays (u, v, w, rho, p, phi, theta) processed from the data file.
    """
    #u = np.zeros((xcells, ycells, zcells, lf))

    file1 = open(fname, 'r')
    print(fname + " file read")
    data = np.loadtxt(fname, skiprows=2)
    k = 0
    for l in range(0, xcells):
        for m in range(0, ycells):
            for n in range(0, zcells):
                u[l, m, n, j] = data[k, 3]
                k += 1

    return u
    
    
def singleRP(case,ord):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))

    folder_case="autotest/caseRPs/"
    folder_out="out/"
    folder_exact="exact/"
    folder_lib="lib"

    fname_config="configure.input"
    #fname_eq="equilibrium.out"
    fname_ini="initial.out"

    folder_out = os.path.join(script_dir, "../"+folder_case+"/"+folder_out)
    folder_exact = os.path.join(script_dir, "../"+folder_case+"/"+folder_exact)
    folder_case = os.path.join(script_dir, "../"+folder_case)
    folder_lib = os.path.join(script_dir, "../"+folder_lib)
    folder_exe = os.path.join(script_dir, "../")


    #Configure the header file for compilation
    backup_file(folder_lib+'/definitions.h')
    modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 2)
    modify_header_file(folder_lib+'/definitions.h', 'ST', 0)
    modify_header_file(folder_lib+'/definitions.h', 'SOLVER', 0)
    modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)
    if case==4: 
        modify_header_file(folder_lib+'/definitions.h', 'MULTICOMPONENT', 1)
        modify_header_file(folder_lib+'/definitions.h', 'MULTI_TYPE', 2)
        
    # COMPILING

    compile_program()
    restore_file(folder_lib+'/definitions.h')

    #Simulation setup
    FinalTime = 0.011
    DumpTime = 0.5
    CFL = 0.4
    Order = ord

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


    #INITIAL CONDITION

    thrs=np.array([1.e-14,0.03,1.2,0.1])

    xc, yc, zc, u, v, w, rho, p, phi, ue, ve, we, rhoe, pe = initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ)
    if case==1: #RP1 Steady
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
        
    if case==4: #RP4
        FinalTime = 0.16
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
        
    if case==3: #RP3
        FinalTime = 0.035
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
        
    if case==5: #RP5 (antes rp3)
        FinalTime = 0.011
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
            
    
      
    # WRITING CONFIGURATION AND INITIAL DATA

    write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, Order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z)
    #write_equilibrium(folder_case, fname_eq, xcells, ycells, zcells, xc, yc, zc, ue, ve, we, rhoe, pe)       
    write_initial(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u, v, w, rho, p, phi)

    # RUNNING

    print("Program is running...")
    print(folder_exe+"./caelum ")
    print(folder_case)
    run_program(folder_exe+"./caelum "+folder_case)

    # READING OUTPUT DATA AND PLOTTING
    files = sorted(glob(folder_out + "/*.out"), key=lambda x: int(re.findall(r'\d+', os.path.basename(x))[0]))
    lf=len(files)
    print(files)

    gamma=1.4
    j=0
    print("Printing figures in folder"+folder_out)
    for fname in files:
        u, v, w, rho, p, phi, theta, E = read_data_euler(fname, xcells, ycells, zcells, lf, gamma, j)   
        j=j+1

    j=-1
    filename = folder_out+"RP_"+str(case)

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
    if case==4: 
        ax[1,1].set_ylabel(r"$1/(1-\gamma)$")

    fig.savefig(filename+".png",dpi=500)

    if case==1:
        L1_error=np.max(u[:,0,0,j])
        
    if 'exactS' in locals():
        interpolated_exact_rho = interp1d(exactS[:, 0], exactS[:, 1], kind='linear', fill_value="extrapolate")
        exact_rho_at_grid = interpolated_exact_rho(xc[:, 0, 0])
        numerical_rho = rho[:, 0, 0, j]
        L1_error = np.mean(np.abs(numerical_rho - exact_rho_at_grid))
    if 'exactSm' in locals():
        interpolated_exact_rho = interp1d(exactSm[:, 0], exactSm[:, 1], kind='linear', fill_value="extrapolate")
        exact_rho_at_grid = interpolated_exact_rho(xc[:, 0, 0])
        numerical_rho = rho[:, 0, 0, j]
        L1_error = np.mean(np.abs(numerical_rho - exact_rho_at_grid))

    print(f"L1 Error: {L1_error}")
    
    if L1_error < thrs[case-1]:
        vf=1
    else:
        vf=0
        
        
    # Remove all .out and .vtk files generated by the simulator
    for f in glob(folder_out + "/*.out") + glob(folder_out + "/*.vtk"):
        os.remove(f)
        
    return vf
    
    
    
    
def caseBubble(ord):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))

    folder_case="autotest/caseBubble/"
    folder_out="out/"
    folder_ref="ref/"
    folder_lib="lib"

    fname_config="configure.input"
    fname_eq="equilibrium.out"
    fname_ini="initial.out"

    folder_out = os.path.join(script_dir, "../"+folder_case+"/"+folder_out)
    folder_ref = os.path.join(script_dir, "../"+folder_case+"/"+folder_ref)
    folder_case = os.path.join(script_dir, "../"+folder_case)
    folder_lib = os.path.join(script_dir, "../"+folder_lib)
    folder_exe = os.path.join(script_dir, "../")


    #Configure the header file for compilation
    backup_file(folder_lib+'/definitions.h')
    modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 2)
    modify_header_file(folder_lib+'/definitions.h', 'ST', 3)
    modify_header_file(folder_lib+'/definitions.h', 'SOLVER', 0)
    modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)

    # COMPILING

    compile_program()
    restore_file(folder_lib+'/definitions.h')

    #Simulation setup
    FinalTime = 800.0
    DumpTime = 100.0
    CFL = 0.45
    Order = ord

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


    #INITIAL CONDITION 

    xc, yc, zc, u, v, w, rho, p, phi, ue, ve, we, rhoe, pe = initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ)
 
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

    # RUNNING

    print("Program is running...")
    run_program(folder_exe+"./caelum "+folder_case)

    # READING OUTPUT DATA AND PLOTTING

    os.remove(folder_out + "list_eq.out")
    files = sorted(glob(folder_out + "/*.out"), key=lambda x: int(re.findall(r'\d+', os.path.basename(x))[0]))
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
        cbar.ax.set_title('θ')
        fig.text(0.15, 0.72, "θ_max="+str(round(np.max(Sth),3)), fontsize=9.5)
        fig.text(0.15, 0.68, "θ_min="+str(round(np.min(Sth),3)), fontsize=9.5)
        
        fig.savefig(filename+".png",dpi=500)
        
        j=j+1
      
    #fname = folder_ref+"/reference.txt"
    #ur, vr, wr, rhor, pr, phir, thetar, Er = read_data_euler(fname, xcells, ycells, zcells, lf, gamma, 0)
    #L1_error = np.mean(np.abs(theta[:,:,:,-1] - thetar[:,:,:,0]))
    
    aux1 = 320.0 - np.max(theta[:,:,:,-1])
    aux2 = np.min(theta[:,:,:,-1]) - 280.0
    
    L1_error = np.max([aux1,aux2])
    
    print(L1_error,aux1,aux2)
    
    if (L1_error < 20) and (L1_error > 1):
        vf=1
    else:
        vf=0
        
    # Remove all .out and .vtk files generated by the simulator
    for f in glob(folder_out + "/*.out") + glob(folder_out + "/*.vtk"):
        os.remove(f)
        
    return vf
    
    
def caseLinear(ord):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))

    folder_case="autotest/caseLinear/"
    folder_out="out/"
    folder_exact="exact/"
    folder_lib="lib"

    fname_config="configure.input"
    #fname_eq="equilibrium.out"
    fname_ini="initial.out"

    folder_out = os.path.join(script_dir, "../"+folder_case+"/"+folder_out)
    folder_exact = os.path.join(script_dir, "../"+folder_case+"/"+folder_exact)
    folder_case = os.path.join(script_dir, "../"+folder_case)
    folder_lib = os.path.join(script_dir, "../"+folder_lib)
    folder_exe = os.path.join(script_dir, "../")


    #Configure the header file for compilation
    backup_file(folder_lib+'/definitions.h')
    modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 0)
    modify_header_file(folder_lib+'/definitions.h', 'ST', 0)
    modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)

    # COMPILING

    compile_program()

    #Simulation setup
    FinalTime = 20.0
    DumpTime = 20.0
    CFL = 0.4
    Order = ord

    #Mesh setup
    xcells = 100
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


    #INITIAL CONDITION

    xc, yc, zc, u, uex, *_  = initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ)
    
    for l in range(0,xcells): 
            for m in range(0,ycells): 
                for n in range(0,zcells):
                    u[l,m,n]=1.0+np.sin(xc[l,m,n]*2.0*math.pi)       
                    uex[l,m,n]=1.0+np.sin((xc[l,m,n]-u_x*FinalTime)*2.0*math.pi) 
      
    # WRITING CONFIGURATION AND INITIAL DATA

    write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, Order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z)      
    write_initial_scalar(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u)

    # RUNNING

    print("Program is running...")
    run_program(folder_exe+"./caelum "+folder_case)

    # READING OUTPUT DATA AND PLOTTING

    files = glob(folder_out+"/*.out")
    lf=len(files)
    print(files)

    j=0
    print("Printing figures in folder"+folder_out)
    u = np.zeros((xcells, ycells, zcells, lf))
    for fname in files:
        u = read_data_scalar(u,fname, xcells, ycells, zcells, lf, j)   
        j=j+1

    j=-1
    filename = folder_out+"linear_scalar_plot"
    
    fig, ax  = plt.subplots(figsize=(8, 6))

    ax.plot(xc[:,0,0],u[:,0,0,j],'o-')
    ax.plot(xc[:,0,0],uex[:,0,0],'-k')
    ax.set_xlabel("x") 
    ax.set_ylabel("u") 
    
    fig.savefig(filename+".png",dpi=500)

    L1_error = np.mean(np.abs(u[:,:,:,-1] - uex[:,:,:]))
    
    print(f"L1 Error: {L1_error}")
    
    restore_file(folder_lib+'/definitions.h')
    
    if L1_error < 0.1*(8-ord):
        vf=1
    else:
        vf=0
        
    # Remove all .out and .vtk files generated by the simulator
    for f in glob(folder_out + "/*.out") + glob(folder_out + "/*.vtk"):
        os.remove(f)
        
    return vf
    
    
def ordersLinear():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    folder_case="autotest/caseLinear/"
    folder_out="out/"
    folder_exact="exact/"
    folder_lib="lib"

    fname_config="configure.input"
    fname_ini="initial.out"

    folder_out = os.path.join(script_dir, "../"+folder_case+"/"+folder_out)
    folder_exact = os.path.join(script_dir, "../"+folder_case+"/"+folder_exact)
    folder_case = os.path.join(script_dir, "../"+folder_case)
    folder_lib = os.path.join(script_dir, "../"+folder_lib)
    folder_exe = os.path.join(script_dir, "../")
    
    backup_file(folder_lib+'/definitions.h')
    modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 0)
    modify_header_file(folder_lib+'/definitions.h', 'ST', 0)
    modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)

    compile_program()
    
    restore_file(folder_lib+'/definitions.h')
      
    orders = [1, 3, 5, 7]
    grid_sizes = [20, 40, 80, 160, 320]
    results = []

    for order in orders:
        for xcells in grid_sizes:
            print(xcells)
            ycells = 1
            zcells = 1
            SizeX = 1.0
            SizeY = 1.0
            SizeZ = 1.0
            
            FinalTime = 5.0
            DumpTime = 20.0
            CFL = 0.01

            Face_1 = 1
            Face_2 = 1
            Face_3 = 1
            Face_4 = 1
            Face_5 = 1
            Face_6 = 1
            
            u_x = 1.0
            u_y = 1.0
            u_z = 1.0


            xc, yc, zc, u, uex, *_  = initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ)
            dx=SizeX/xcells
            for l in range(xcells):  
                for m in range(ycells): 
                    for n in range(zcells):
                        x1 = xc[l,m,n] + dx/2.0
                        x2 = xc[l,m,n] - dx/2.0
                        u[l,m,n] =  1.0 + 0.5*(-np.cos(2.0*math.pi*x1) + np.cos(2.0*math.pi*x2))/(dx*2.0*math.pi)      
                        uex[l,m,n] =  1.0 + 0.5*(-np.cos(2.0*math.pi*(x1-u_x*FinalTime)) + np.cos(2.0*math.pi*(x2-u_x*FinalTime)))/(dx*2.0*math.pi)   

            write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z)      
            write_initial_scalar(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u)

            print("Program is running...")
            run_program(folder_exe + "./caelum " + folder_case)

            files = sorted(glob(folder_out + "/*.out"), key=lambda x: int(re.findall(r'\d+', os.path.basename(x))[0]))
            lf = len(files)
            print(files)

            j = 0
            print("Printing figures in folder " + folder_out)
            u = np.zeros((xcells, ycells, zcells, lf))
            for fname in files:
                u = read_data_scalar(u,fname, xcells, ycells, zcells, lf, j)   
                j += 1

            L1_error = np.sum(np.abs(u[:,:,:,-1] - uex[:,:,:]))*dx
            #print(f"Order {order}, Grid {xcells}: L1 Error = {L1_error}")
            results.append((order, xcells, L1_error))

            
    
    # Compute the orders of accuracy
    orders_accuracy = {}
    mean_order = np.zeros(len(orders))
    ii=0
    for order in orders:
        L1_errors = [result[2] for result in results if result[0] == order]
        cell_sizes = [1.0 / size for size in grid_sizes]
        order_accuracy = []
        for jj in range(len(cell_sizes) - 1):
            order_acc = np.log(L1_errors[jj] / L1_errors[jj + 1]) / np.log(cell_sizes[jj] / cell_sizes[jj + 1])
            order_accuracy.append(order_acc)
            mean_order[ii] += order_acc/(len(cell_sizes)-1)

        orders_accuracy[order] = order_accuracy
        ii+=1
    
    
    # Print results and orders of accuracy to a text file
    output_file_path = os.path.join(folder_out, "results.txt")
    with open(output_file_path, "w") as f:
        for order, xcells, error in results:
            f.write(f"Order {order}, Grid {xcells} cells: L1 Error = {error}\n")
        for order in orders_accuracy:
            f.write(f"Order {order} accuracy: {orders_accuracy[order]}\n")
    
    order_errors  = {order: [] for order in orders}
    grid_errors   = {grid : [] for grid  in grid_sizes}
    for order in orders:
        for xcells in grid_sizes:
            for result in results:
                if result[0] == order and result[1] == xcells:
                    order_errors[order].append(result[2])
                    grid_errors[xcells].append(result[2])

    # Plotting the data
    fig, (ax1,ax2)  = plt.subplots(1,2,figsize=(12, 4))
    for order in orders:
        ax1.plot(grid_sizes, order_errors[order], 'o-', label=f'Order {order}')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel("Number of cells") 
    ax1.set_ylabel("$L_1$") 
    ax1.legend()
    for xcells in grid_sizes:
        ax2.plot(orders, grid_errors[xcells], 'o-',label=f'{xcells} cells')
    ax2.set_yscale('log')
    ax2.set_xlabel("Order") 
    ax2.set_ylabel("$L_1$") 
    ax2.legend()
    
    x0 = 20
    y0 = 1e-9
    x_range = np.array([40, 200])
    y0_vals = np.array([4.e-1, 4.e-2, 1.e-3, 1.e-4])
    theoretical_orders = [1, 3, 5, 7]
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    #ax1.plot(x_range, [y0,y0], '-', lw=0.5, color="k")    
    #ax1.plot([60,60], [y0,y0 * (60 / x0) ** (-7)], '-', lw=0.5, color="k")
    jj=0
    for p, color in zip(theoretical_orders, colors):
        y_range = y0_vals[jj] * (x_range / x0) ** (-p)
        ax1.plot(x_range, y_range, '--', lw=0.75, color=color)
        jj+=1
    
    
    filename = folder_out+"convergences"
    fig.savefig(filename+".png",dpi=500)
    
    print("\nThe mean orders are: ", mean_order) 
    comp = np.abs(mean_order - np.array(orders))
    mean_ord_error = np.mean(comp)
    print("The mean difference in orders is: ", mean_ord_error, "\n")
    
    if mean_ord_error < 1.0:
        vf=1
    else:
        vf=0
        
    # Remove all .out and .vtk files generated by the simulator
    for f in glob(folder_out + "/*.out") + glob(folder_out + "/*.vtk"):
        os.remove(f)

    return vf
    
    
def ordersEuler():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    folder_case="autotest/caseAccuracy/"
    folder_out="out/"
    folder_exact="exact/"
    folder_lib="lib"

    fname_config="configure.input"
    fname_ini="initial.out"

    folder_out = os.path.join(script_dir, "../"+folder_case+"/"+folder_out)
    folder_exact = os.path.join(script_dir, "../"+folder_case+"/"+folder_exact)
    folder_case = os.path.join(script_dir, "../"+folder_case)
    folder_lib = os.path.join(script_dir, "../"+folder_lib)
    folder_exe = os.path.join(script_dir, "../")

    orders = [1, 3, 5, 7]
    grid_sizes = [20, 40, 80, 160, 320]
    results = []

    for order in orders:
        for xcells in grid_sizes:
            print(xcells)
            ycells = 1
            zcells = 1
            SizeX = 1.0
            SizeY = 1.0
            SizeZ = 1.0
            
            FinalTime = 5.0
            DumpTime = 20.0
            CFL = 0.01

            Face_1 = 1
            Face_2 = 1
            Face_3 = 1
            Face_4 = 1
            Face_5 = 1
            Face_6 = 1
            
            u_x = 1.0
            u_y = 1.0
            u_z = 1.0

            backup_file(folder_lib+'/definitions.h')
            modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 2)
            modify_header_file(folder_lib+'/definitions.h', 'ST', 0)
            modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)

            compile_program()
            
            uex = np.zeros((xcells, ycells, zcells))
            xc, yc, zc, u, v, w, rho, p, phi, ue, ve, we, rhoe, pe = initialize_variables(xcells, ycells, zcells, SizeX, SizeY, SizeZ)
            dx=SizeX/xcells
            for l in range(xcells):  
                for m in range(ycells): 
                    for n in range(zcells):
                        x1 = xc[l,m,n] + dx/2.0
                        x2 = xc[l,m,n] - dx/2.0
                        rho[l,m,n] =  1.0 + 0.5*(-np.cos(2.0*math.pi*x1) + np.cos(2.0*math.pi*x2))/(dx*2.0*math.pi)
                        p  [l,m,n]=1.0
                        u  [l,m,n]=u_x
                        v  [l,m,n]=0.0
                        w  [l,m,n]=0.0
                        phi[l,m,n]=1.0                        
                        uex[l,m,n] =  1.0 + 0.5*(-np.cos(2.0*math.pi*(x1-u_x*FinalTime)) + np.cos(2.0*math.pi*(x2-u_x*FinalTime)))/(dx*2.0*math.pi)   

            write_config(folder_case, fname_config, FinalTime, DumpTime, CFL, order, xcells, ycells, zcells, SizeX, SizeY, SizeZ, Face_1, Face_2, Face_3, Face_4, Face_5, Face_6, u_x, u_y, u_z)      
            write_initial(folder_case, fname_ini, xcells, ycells, zcells, xc, yc, zc, u, v, w, rho, p, phi)

            print("Program is running...")
            run_program(folder_exe + "./caelum " + folder_case)

            files = glob(folder_out + "/*.out")
            lf = len(files)
            print(files)

            gamma=1.4
            j = 0
            print("Printing figures in folder " + folder_out)
            for fname in files:
                u, v, w, rho, p, phi, theta, E = read_data_euler(fname, xcells, ycells, zcells, lf, gamma, j)  
                j += 1

            #print(rho[:,0,0,-1])
            #print(uex[:,0,0])
            error=np.abs(rho[:,0,0,-1] - uex[:,0,0])
            #print(error)
            L1_error = np.sum(error)*dx
            print(f"Order {order}, Grid {xcells}: L1 Error = {L1_error}")
            results.append((order, xcells, L1_error))

            restore_file(folder_lib+'/definitions.h')
    
    # Compute the orders of accuracy
    orders_accuracy = {}
    for order in orders:
        L1_errors = [result[2] for result in results if result[0] == order]
        cell_sizes = [1.0 / size for size in grid_sizes]
        order_accuracy = []
        for jj in range(len(cell_sizes) - 1):
            order_acc = np.log(L1_errors[jj] / L1_errors[jj + 1]) / np.log(cell_sizes[jj] / cell_sizes[jj + 1])
            order_accuracy.append(order_acc)
        orders_accuracy[order] = order_accuracy
    
    for order in orders_accuracy:
        print(f"Order {order} accuracy: {orders_accuracy[order]}")

    return results, orders_accuracy
    