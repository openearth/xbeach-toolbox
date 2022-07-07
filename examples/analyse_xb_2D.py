import matplotlib.pyplot as plt
from scripts.xb_analyse import XBeachModelAnalysis

r03 = XBeachModelAnalysis('r03', r'c:\Users\marliesvanderl\phd\modeling\rvw_phzd\runs\r03')

# load the xbeach model set-up
r03.load_model_setup()

# get an idea of the grid shape
r03.load_modeloutput('H')
ny, nx = r03.var['globalx'].shape

# check the boundary conditions for the tide
fig, ax = r03.fig_check_tide_bc()

# change coordinates of plots to local coordinates:
r03.set_plot_localcoords(True)

# only plot a certain Area Of Interest of the complete grid
# r03.set_aoi([20,445,20,220])

# example usage map plotting
fig, ax = r03.fig_map_diffvar('zb', '$\Delta z_b$ [m]', itend=9, it0=0)
fig, ax = r03.fig_map_var('H','Hm0 [m]')

# r03.load_modeloutput('point_H')

# load point model data for comparison with field data
data = r03.get_modeloutput_by_station('point_H','xylem')
plt.figure()
plt.plot(data[:,0],data[:,1], label = 'xb xylem')
