# Notes on xBeach-ToolBox

## Folders

### General

1) deg2uv.py

This script contains one function that converts degrees to u,v vectors and the reverse

Functions:
    * deg2uv - converts degrees to u, v vectors
    * uv2deg - converts u, v vectors to degrees

2) executing_runs.py
   * xb_run_script_win     - Create batch script to run simulations. (Not sure exactly what this is doing. Especially what the N parameter does)
   * generate_batch_script - generates a batch script for a xbeach model
   * run_batch_script      - Runs a batch script. current not working all the way
  
3) geometry.py
   * rotate                - rotates a vecotor through the origin that is set by the user
   * rotate_grid           - Rotates a grid
   * grid_world2local      - Figures out the rotation of a grid and then returns it in local coordinates
   * samples_to_local_grid - shifts a vector and then rotates it
   * in_polygon            - Checks whether a set of coordinates are in a polygon
   * get_points_in_polygon - Returns True for points in polygon and Returns False for points outside of the polygon - May replace in_polygon in the future
   * path_distance         - Computes the distance along a polyline
   * 

4) visualize_mesh.py
   * plot_mesh                - Plot the mesh based on the meshgrid regular grid - **Uses **kwargs to feed inputs into the plotting functions. This is a good idea**
   * write_mesh_to_shp        - Convert mesh to shapefile
   * grid_parameter_to_shape  - Saves a paramter defined at a grid point to a shapefile
  
5) wave_functions
   * dispersion                - Computes the wave nuber given a radial frequency and water depth
   * wavecelerity              - Computes the group velocity and wave celerity ratio
   * directional_spread_kuik   - Determine mean wave direction and directional spreading
   * celertiry_ratio_equals_09 - Function to find water depth for which the n ratio equal 0.9
   * offshore_depth            - Compute the required offshore depth to correctly force the waves

### grid

6) creation.py
   * xgrid                 - Compute spatially varying grid based on the local wave length
   * ygrid                 - function to setup a basic ygrid array
   * grid_transition       - Function for grid transition from one cell to the other. (Still not sure what this means)

7) extension.py
   * lateral_extend        - Extend the model domain at both laterial sides with n number of cells
   * seaward_extend        - Compute the seaward extend of the bathymery based on an artificial  slope and required offshore depth

8) grid_check.py
   * diff                  - Computes the difference in subsequent items
   * smoothness            - Perform grid checks based on D3D-quickin procedure for smoothness in m and n-direction and stores figures
   * aspect_ratio          - Perform grid checks based on D3D-quickin procedure
   * orthogonality         - Perform grid checks based on D3D-quickin procedure for orthogonality and stores figures

9) refinement.py
   * grid_refine_grid      - refines the grid with the factor xfactor in xdirection and yfactor in y direction