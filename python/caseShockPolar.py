import os
import math                    
import numpy as np             
import matplotlib.pyplot as plt 
from matplotlib.patches import Polygon
from glob import glob
import re
import imageio

from utils import (
    modify_header_file, write_config, write_initial,
    write_equilibrium, backup_file, restore_file,
    compile_program, run_program, initialize_variables,
    read_data_euler, write_solid_cells
)

# =========================================================
# SHOCK POLAR ANALYTICAL RELATION
# =========================================================

def theta_from_beta(beta, M, gamma=1.4):
    num = 2 / np.tan(beta) * (M**2 * np.sin(beta)**2 - 1)
    den = M**2 * (gamma + np.cos(2*beta)) + 2
    return np.arctan(num / den)

# =========================================================
# DOMAIN SETUP
# =========================================================

folder_case = "caseShockPolar/"
script_dir = os.path.dirname(os.path.abspath(__file__))

folder_out = os.path.join(script_dir, "../run/"+folder_case+"/out/")
folder_case = os.path.join(script_dir, "../run/"+folder_case)
folder_lib = os.path.join(script_dir, "../lib")
folder_exe = os.path.join(script_dir, "../")

os.makedirs(folder_case, exist_ok=True)
os.makedirs(folder_out, exist_ok=True)

# =========================================================
# PARAMETER SWEEP (Mach + theta)
# =========================================================

Mach_list = [1.5, 2.0]
theta_list_deg = np.linspace(15, 20, 2)

results = {M: {"theta": [], "beta": []} for M in Mach_list}

gamma = 1.4

# =========================================================
# BACKUP COMPILATION FILE
# =========================================================

backup_file(folder_lib+'/definitions.h')

modify_header_file(folder_lib+'/definitions.h', 'NTHREADS', 48)
modify_header_file(folder_lib+'/definitions.h', 'TYPE_REC', 0)
modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 2)
modify_header_file(folder_lib+'/definitions.h', 'ST', 0)
modify_header_file(folder_lib+'/definitions.h', 'SOLVER', 0)
modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)

compile_program()
restore_file(folder_lib+'/definitions.h')

# =========================================================
# LOOP OVER MACH AND THETA
# =========================================================

for Ms in Mach_list:
    for theta_deg in theta_list_deg:

        print(f"\n==============================")
        print(f"Running M={Ms}, theta={theta_deg}")
        print(f"==============================")

        # -------------------------
        # clean output
        # -------------------------
        for f in glob(folder_out + "/*.out") + glob(folder_out + "/*.vtk"):
            os.remove(f)

        # -------------------------
        # geometry
        # -------------------------
        theta = np.deg2rad(theta_deg)

        xcells = 100
        ycells = 250
        zcells = 1

        SizeX = 0.1
        SizeY = 0.25
        SizeZ = SizeY/ycells

        xc, yc, zc, u, v, w, rho, p, phi, ue, ve, we, rhoe, pe = initialize_variables(
            xcells, ycells, zcells, SizeX, SizeY, SizeZ
        )

        sld = np.zeros((xcells, ycells, zcells))

        # -------------------------
        # physical state
        # -------------------------
        pL = 0.05e6
        TL = 285.0
        Rgas = 287.0

        rhoL = pL/(Rgas*TL)
        aL = np.sqrt(gamma*pL/rhoL)

        uL = Ms * aL
        vL = 0.0

        rhoR = rhoL
        pR = pL
        uR = 0.0
        vR = 0.0

        # -------------------------
        # geometry triangle
        # -------------------------
        xa = 0.02
        ya = SizeY/2
        length_x = 0.1

        slope = np.tan(theta)

        x_base = xa + length_x

        x0, y0 = xa, ya
        x1, y1 = x_base, ya + slope*length_x
        x2, y2 = x_base, ya - slope*length_x

        v0 = np.array([x0, y0])
        v1 = np.array([x1, y1])
        v2 = np.array([x2, y2])

        xs = 0.01

        # -------------------------
        # INITIAL CONDITION
        # -------------------------
        for l in range(xcells):
            for m in range(ycells):
                for n in range(zcells):

                    x = xc[l,m,n]
                    y = yc[l,m,n]
                    p_point = np.array([x,y])

                    if x < xs:
                        p[l,m,n] = pL
                        rho[l,m,n] = rhoL
                        u[l,m,n] = uL
                        v[l,m,n] = vL
                    else:
                        p[l,m,n] = pR
                        rho[l,m,n] = rhoR
                        u[l,m,n] = uR
                        v[l,m,n] = vR

                    # signed distance (optional)
                    sld[l,m,n] = 0.0

        # -------------------------
        # WRITE INPUT FILES
        # -------------------------
        write_config(folder_case, "configure.input",
                     0.002, 0.0004, 0.4, 3,
                     xcells, ycells, zcells,
                     SizeX, SizeY, SizeZ,
                     3,3,3,3,1,1,
                     1,1,1)

        write_initial(folder_case, "initial.out",
                     xcells, ycells, zcells,
                     xc, yc, zc,
                     u, v, w, rho, p, phi)

        write_solid_cells(folder_case, "solid_cells.input",
                         xcells, ycells, zcells,
                         xc, yc, zc, sld)

        # -------------------------
        # RUN SIMULATION
        # -------------------------
        run_program(folder_exe+"./caelum "+folder_case)

        # =====================================================
        # POSTPROCESS: extract beta
        # =====================================================

        files = sorted(glob(folder_out + "/*.out"),
                       key=lambda x: int(re.findall(r'\d+', x)[0]))

        for fname in files[-1:]:  # last timestep only

            u, v, w, rho, p, phi, theta_field, E = read_data_euler(
                fname, xcells, ycells, zcells, 1, gamma, 0
            )

            xp = xc[:,0,0]
            yp = yc[0,:,0]

            Srho = np.transpose(rho[:,:,0,0])

            dx = xp[1]-xp[0]
            dy = yp[1]-yp[0]

            drho_dy, drho_dx = np.gradient(Srho, dy, dx)
            grad = np.sqrt(drho_dx**2 + drho_dy**2)

            threshold = 0.1*np.max(grad)

            x_shock = []
            y_shock = []

            for j in range(grad.shape[0]):
                row = grad[j,:]
                valid = row > threshold

                if np.any(valid):
                    idx = np.argmax(row*valid)
                    x_shock.append(xp[idx])
                    y_shock.append(yp[j])

            if len(x_shock) > 5:

                coeff = np.polyfit(x_shock, y_shock, 1)
                beta = np.arctan(coeff[0])
                beta_deg = np.rad2deg(beta)

                results[Ms]["theta"].append(theta_deg)
                results[Ms]["beta"].append(beta_deg)

                print(f"β = {beta_deg:.2f}")

# =========================================================
# FINAL SHOCK POLAR PLOT
# =========================================================

plt.figure(figsize=(8,6))

beta_range = np.linspace(np.deg2rad(1), np.deg2rad(89), 200)

for M in Mach_list:
    theta_vals = theta_from_beta(beta_range, M, gamma)

    valid = theta_vals > 0

    plt.plot(np.rad2deg(theta_vals[valid]),
             np.rad2deg(beta_range[valid]),
             label=f"Analítico M={M}")

for M in Mach_list:
    plt.scatter(results[M]["theta"],
                results[M]["beta"],
                label=f"CFD M={M}")

plt.xlabel("θ (deg)")
plt.ylabel("β (deg)")
plt.title("Shock Polar θ–β–M")
plt.grid()
plt.legend()
plt.tight_layout()

plt.savefig(folder_out + "shock_polar.png", dpi=300)
plt.show()