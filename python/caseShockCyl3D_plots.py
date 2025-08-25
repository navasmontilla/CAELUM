import os
from glob import glob
import re
import pyvista as pv
import imageio
import numpy as np


folder_case="caseShockCyl3D/"
script_dir = os.path.dirname(os.path.abspath(__file__))
folder_out="out/"
folder_out = os.path.join(script_dir, "../run/"+folder_case+"/"+folder_out)

# Start the virtual framebuffer (Xvfb) to enable off-screen rendering
pv.start_xvfb()

# Find output files
vtk_files = glob(folder_out + "/*.vtk")

# Filter files that contain digits in their filenames
vtk_files_with_digits = [f for f in vtk_files if re.search(r'\d+', os.path.basename(f))]

# Sort the filtered files by the first sequence of digits found in the filenames
files = sorted(vtk_files_with_digits, key=lambda x: int(re.findall(r'\d+', os.path.basename(x))[0]))
lf=len(files)
print(files)

j=0
images = [] 
print("Printing figures in folder"+folder_out)

for fname in files:
    
    vtk_file=fname
    data = pv.read(vtk_file)
    grid = data.compute_cell_sizes().cell_data_to_point_data()
    plotter = pv.Plotter(off_screen=True,window_size=[2200, 1400])

    #opacity = [0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.90, 1.0]  
    #plotter.add_volume(grid, scalars="rho", cmap="RdGy_r",opacity=opacity,opacity_unit_distance=0.1,clim=[0.2, 2.0]) #antes 0.1
    
    #contour_values = [0.5,1.0]
    #contours = grid.contour(isosurfaces=contour_values, scalars="rho") 
    #if contours.n_points > 0:
    #    plotter.add_mesh(contours, cmap="coolwarm", scalars="rho", opacity=0.4)

    # Compute the gradient of the rho field
    gradient = grid.compute_derivative(scalars="rho", gradient=True)
    
    # Manually compute the magnitude of the gradient vector field
    gradient_array = gradient.point_data['gradient']
    gradient_magnitude = np.linalg.norm(gradient_array, axis=1)
    
    # Add the gradient magnitude as a scalar field to the grid
    grid.point_data['gradient_magnitude'] = gradient_magnitude

    # Add the Schlieren effect visualization
    #plotter.add_mesh(grid, scalars="gradient_magnitude", cmap="bone_r", opacity=0.6)
    
    opacity = [0.0, 0.9, 0.99, 0.99, 1.0]  
    plotter.add_volume(grid, scalars="gradient_magnitude", cmap="bone_r",opacity=opacity,opacity_unit_distance=0.02) #antes 0.1

    #contour_values = [100,101]
    #contours = grid.contour(isosurfaces=contour_values, scalars="gradient_magnitude") 
    #if contours.n_points > 0:
    #    plotter.add_mesh(contours, color="tab:red", opacity=0.4)
    
    plotter.camera_position = [
        (3.2,-1.5, 1.5),  # Camera position (x, y, z)
        (0.75, 0.5, 0.25),  # Focal point (center of the object)
        (0, 0, 1),  # View up vector (defines the up direction)
    ]

    plotter.show_axes()
    plotter.add_bounding_box(color='black') 
    
    image_path = fname+".png"
    plotter.show(screenshot=image_path)
    img = plotter.screenshot(return_img=True)
    images.append(img)
    
    
    j=j+1
         

gif_path = os.path.join(folder_out, "animation.gif")
imageio.mimsave(gif_path, images, fps=1, loop=0)  # Adjust the duration as needed

