from scripts.xb_analyse import XBeachModelAnalysis

r03 = XBeachModelAnalysis('r03', r'p:\11208248-xb-verkenning-hrd\modeling\r18')

r03.load_model_setup()

r03.read_modeloutput('H')

fig, ax = r03.fig_check_tide_bc()
