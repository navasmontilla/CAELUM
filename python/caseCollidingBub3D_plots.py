import os
from glob import glob
import re
import pyvista as pv
import imageio


folder_case="caseCollidingBub3D/"
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

    opacity = [1.0, 0.9, 0.0, 0.0, 0.0, 0.9, 1.0]  
    plotter.add_volume(grid, scalars="theta", cmap="RdBu_r",opacity=opacity,opacity_unit_distance=200,clim=[290, 310]) 

    plotter.camera_position = [
        (30000, 30000, 12000),  # Camera position (x, y, z)
        (10000, 10000, 5000),  # Focal point (center of the object)
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
imageio.mimsave(gif_path, images, duration=8, loop=0)  # Adjust the duration as needed

