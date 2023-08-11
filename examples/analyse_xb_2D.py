import matplotlib.pyplot as plt
import sys
import os

## import xbeach tools
sys.path.append(os.path.abspath(os.path.join('..')))

from xbTools import XBeachModelAnalysis



results = XBeachModelAnalysis('r03', 'p:\\xbeach\skillbed\\test_menno_v2\\run\Saint_Trojan\\')


## additional settings (not needed to specify)
results.set_save_fig(False)
# plot in local coords instead of global coords
results.set_plot_localcoords(False)
# plot in kilometers instead of m
results.set_plot_km_coords(False)
# set starttime, if not specified time from netcdf is used
#results.set_globalstarttime('2021-10-11T13:00:00')
# change units
results.set_unitdict({'zb':['m+msl']}) 


# load the xbeach model set-up
results.load_model_setup()




zs = results.load_modeloutput('zs')

ny, nx = results.var['globalx'].shape

results.load_modeloutput('point_zs')


zs = results.get_modeloutput('zs')
point_zs = results.get_modeloutput('point_zs')

[t, zs] = results.get_modeloutput_by_station('zs','point001')


results.fig_check_tide_bc()



## change coordinates of plots to local coordinates:
#results.set_plot_localcoords(True)
## only plot a certain Area Of Interest of the complete grid
#results.set_aoi([20,445,20,220])

# example usage map plotting
fig, ax = results.fig_map_var('u','m/s')


fig, ax = results.fig_map_diffvar('zb', '$\Delta z_b$ [m]', itend=9, it0=0)

fig, ax = results.fig_cross_var('zs', 5, iy=5, coord=None, plot_ref_bathy=True, zmin=-25)

fig, ax = results.fig_profile_change(iy=5)




