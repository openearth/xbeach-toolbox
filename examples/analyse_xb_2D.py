from scripts.xb_analyse import XBeachModelAnalysis

r03 = XBeachModelAnalysis('r03', r'c:\Users\marliesvanderl\phd\modeling\rvw_phzd\runs\r01')

r03.load_model_setup()

r03.read_modeloutput('H')
ny, nx = r03.var['globalx'].shape

fig, ax = r03.fig_check_tide_bc()

r03.set_plot_localcoords(True)
r03.set_aoi([20,445,20,220])

fig, ax = r03.fig_map_diffvar('zb', '$\Delta z_b$ [m]', itend=20, it0=0)
fig, ax = r03.fig_map_var('H','Hm0 [m]')