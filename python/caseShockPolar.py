import os
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import re

from utils import (
    modify_header_file, write_config, write_initial,
    backup_file, restore_file, compile_program,
    run_program, initialize_variables,
    read_data_euler, write_solid_cells
)

# =========================================================
# SIGNED DISTANCE TRIANGLE
# =========================================================
def udTriangle2D(p, a, b, c):

    def seg_sdf(pv, v0, v1):
        v = v1 - v0
        w = pv - v0
        t = np.clip(np.dot(w, v)/np.dot(v,v), 0.0, 1.0)
        proj = v0 + t*v
        return np.linalg.norm(pv - proj)

    d = min(seg_sdf(p, a, b),
            seg_sdf(p, b, c),
            seg_sdf(p, c, a))

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

    return -d if (u >= 0 and v >= 0 and u+v <= 1) else d


# =========================================================
# SHOCK POLAR ANALYTICAL
# =========================================================
def theta_from_beta(beta, M, gamma=1.4):
    num = 2/np.tan(beta)*(M**2*np.sin(beta)**2 - 1)
    den = M**2*(gamma + np.cos(2*beta)) + 2
    return np.arctan(num/den)


# =========================================================
# PATHS
# =========================================================
folder_case_name = "caseShockPolar/"
script_dir = os.path.dirname(os.path.abspath(__file__))

folder_case = os.path.join(script_dir, "../run/"+folder_case_name)
folder_out  = os.path.join(folder_case, "out/")
folder_lib  = os.path.join(script_dir, "../lib")
folder_exe  = os.path.join(script_dir, "../")

os.makedirs(folder_case, exist_ok=True)
os.makedirs(folder_out, exist_ok=True)


# =========================================================
# SIMULATION PARAMETERS 
# =========================================================

# --- Time ---
FinalTime = 0.002
DumpTime  = 0.0004
CFL       = 0.4
Order     = 3

# --- Mesh ---
xcells = 100
ycells = 250
zcells = 1

SizeX = 0.1
SizeY = 0.25
SizeZ = SizeY / ycells

# --- Boundary conditions ---
Face_1 = 3
Face_2 = 3
Face_3 = 3
Face_4 = 3
Face_5 = 1
Face_6 = 1

# --- Transport ---
u_x = 1.0
u_y = 1.0
u_z = 1.0


# =========================================================
# PARAM SWEEP
# =========================================================
Mach_list = [2.0,2.5]
#theta_list_deg = np.linspace(12, 20, 2)
theta_dict = {
    2.0: np.linspace(12, 22, 4),
    2.5: np.linspace(12, 29, 5)
}

results = {M: {"theta": [], "beta": []} for M in Mach_list}

gamma = 1.4


# =========================================================
# COMPILATION
# =========================================================
backup_file(folder_lib+'/definitions.h')

modify_header_file(folder_lib+'/definitions.h', 'NTHREADS', 48)
modify_header_file(folder_lib+'/definitions.h', 'EQUATION_SYSTEM', 2)
modify_header_file(folder_lib+'/definitions.h', 'READ_INITIAL', 1)
modify_header_file(folder_lib+'/definitions.h', 'ALLOW_SOLIDS', 3)

compile_program()
restore_file(folder_lib+'/definitions.h')


# =========================================================
# LOOP
# =========================================================
for Ms in Mach_list:
    theta_list_deg = theta_dict[Ms]
    for theta_deg in theta_list_deg:

        print(f"\n=== M={Ms}, θ={theta_deg} ===")

        # limpiar salida
        for f in glob(folder_out + "/*.out"):
            os.remove(f)

        # -------------------------
        # INIT ARRAYS
        # -------------------------
        xc, yc, zc, u, v, w, rho, p, phi, *_ = initialize_variables(
            xcells, ycells, zcells, SizeX, SizeY, SizeZ
        )

        sld = np.zeros((xcells, ycells, zcells))

        # -------------------------
        # FLOW STATE
        # -------------------------
        pL = 0.05e6
        TL = 285.0
        Rgas = 287.0

        rhoL = pL/(Rgas*TL)
        aL = np.sqrt(gamma*pL/rhoL)

        uL = Ms * aL
        vL = 0.0

        # -------------------------
        # GEOMETRY (CUÑA)
        # -------------------------
        theta = np.deg2rad(theta_deg)

        xa = 0.02
        ya = SizeY/2
        L  = 0.06

        slope = np.tan(theta)

        x0, y0 = xa, ya
        x1, y1 = xa + L, ya + slope*L
        x2, y2 = xa + L, ya - slope*L

        v0 = np.array([x0, y0])
        v1 = np.array([x1, y1])
        v2 = np.array([x2, y2])

        # -------------------------
        # INITIAL CONDITION
        # -------------------------
        for l in range(xcells):
            for m in range(ycells):
                for n in range(zcells):

                    x = xc[l,m,n]
                    y = yc[l,m,n]

                    rho[l,m,n] = rhoL
                    p[l,m,n]   = pL
                    u[l,m,n]   = uL
                    v[l,m,n]   = vL

                    sld[l,m,n] = udTriangle2D(
                        np.array([x,y]), v0, v1, v2
                    )

        # -------------------------
        # WRITE FILES
        # -------------------------
        write_config(folder_case, "configure.input",
                     FinalTime, DumpTime, CFL, Order,
                     xcells, ycells, zcells,
                     SizeX, SizeY, SizeZ,
                     Face_1, Face_2, Face_3, Face_4, Face_5, Face_6,
                     u_x, u_y, u_z)

        write_initial(folder_case, "initial.out",
                      xcells, ycells, zcells,
                      xc, yc, zc,
                      u, v, w, rho, p, phi)

        write_solid_cells(folder_case, "solid_cells.input",
                          xcells, ycells, zcells,
                          xc, yc, zc, sld)

        # -------------------------
        # RUN
        # -------------------------
        run_program(folder_exe+"./caelum "+folder_case)

        # =====================================================
        # POSTPROCESS
        # =====================================================
        files = sorted(glob(folder_out + "/*.out"),
                       key=lambda x: int(re.findall(r'\d+', x)[0]))

        fname = files[-1]

        u, v, w, rho, p, phi, _, _ = read_data_euler(
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

        x_shock = np.array(x_shock)
        y_shock = np.array(y_shock)

        # filtro espacial
        mask = (
            (x_shock > x0) &
            (y_shock > SizeY*0.65) &
            (y_shock < SizeY*0.85)
        )

        x_shock = x_shock[mask]
        y_shock = y_shock[mask]

        if len(x_shock) > 10:

            coeff = np.polyfit(x_shock, y_shock, 1)
            beta = np.arctan(coeff[0])
            beta_deg = np.rad2deg(beta)

            results[Ms]["theta"].append(theta_deg)
            results[Ms]["beta"].append(beta_deg)

            print(f"β_num = {beta_deg:.2f}")


# =========================================================
# FINAL PLOT
# =========================================================
plt.figure(figsize=(7,5))

beta_range = np.linspace(np.deg2rad(1), np.deg2rad(89), 400)

for M in Mach_list:
    theta_vals = theta_from_beta(beta_range, M, gamma)
    valid = theta_vals > 0

    plt.plot(np.rad2deg(theta_vals[valid]),
             np.rad2deg(beta_range[valid]),
              label=f"Analytical M={M}")

for M in Mach_list:
    plt.scatter(results[M]["theta"],
                results[M]["beta"],
                label=f"Numerical (IBM) M={M}")

plt.xlabel("θ (deg)")
plt.ylabel("β (deg)")
plt.title("Shock Polar θ–β–M")
plt.grid()
plt.legend()
plt.tight_layout()

plt.savefig(folder_out + "shock_polar.png", dpi=300)
plt.show()