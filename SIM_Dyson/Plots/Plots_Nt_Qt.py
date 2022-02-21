import logging
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.ticker import FormatStrFormatter
pl.style.use('classic')

import Functions as Fun

pl.rcParams["figure.figsize"] = [15, 11]

args = Fun.argument_parsing()
loglevel = args.log
if not args.silent :
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(format='%(asctime)s : %(message)s', level=numeric_level)

N = args.N
N_iter = args.Niter
Nsims = args.Nsim
dt = args.dt
lf = args.lf
T = N_iter*dt
Qmax = N
beta_min = args.beta_min
beta_max = args.beta_max
delta_beta = args.delta_beta

logging.info("Plot program called with arguments:\n %s", args)

betas = np.arange(beta_min, beta_max+delta_beta, delta_beta)
time = Fun.time_df(dt, N_iter, steps = 100)

logging.info("Starting to make the plots for N(t)")

Tit = r'$N_0 = $' + str(N) + r', $l_f = $' + str(lf)
figname = 'Ninit' + str(N) + '_lf' + str(lf) + '_dt' + str(dt) + '.pdf' # To save plot

N_normalize = mcolors.Normalize(vmin = 0.1, vmax = 1.0)
N_colormap = cm.get_cmap('jet', len(betas)+1)

N_avg_t = pl.figure()

ax = pl.axes()
ax = pl.axes(xscale = 'log', yscale = 'log')
ax.grid()
ax.set_xlim(time['t'].min(), time['t'].max())
pl.xlabel(r'$\ln(\delta t/f)$', fontsize = 23)
pl.ylabel(r'$\ln\left(\langle N(t) \rangle\right)$', fontsize = 23)
pl.title(Tit, fontsize = 25, pad = 15)
pl.xticks(fontsize = 16)
pl.yticks(fontsize = 16)
locs, labels = pl.yticks()
pl.yticks(np.arange(1, N + 2, 2), np.arange(1, N + 2, 2))
for beta in reversed(betas):
	filenamebase = Fun.build_filenamebase(args, beta)
	logging.debug("filenamebase = %s", filenamebase)

	Ndf_filename = './' + filenamebase+'/Ndt.csv'
	logging.info("Starting read data from file %s", Ndf_filename)

	N_df = pd.read_csv(Ndf_filename)

	lab = r'$\beta = $' + str(beta)
	ax.set_ylim(N_df['N_avg'].min(), N_df['N_avg'].max() + 2)
	#ax.errorbar(N_df['t'], N_df['N_avg'], yerr = N_df['std_N'], label = lab, linewidth = 0.75, ms = 1.5)
	pl.plot(N_df['t'], N_df['N_avg'], label = lab, color = N_colormap(N_normalize(beta)), linewidth = 0.9)
	pl.fill_between(N_df['t'], N_df['N_avg'] + N_df['std_N'], N_df['N_avg'] - N_df['std_N'],#
		color = N_colormap(N_normalize(beta)), alpha = 0.2, linewidth = 0.75)

pl.subplots_adjust(left = 0.095, bottom = 0.075, right = 0.95, top = 0.925)
pl.legend(loc = 'lower left', fontsize = 18)
N_avg_t.savefig(figname, dpi = 500)
pl.close(N_avg_t)


for beta in betas:
	filenamebase = Fun.build_filenamebase(args, beta)
	logging.debug("filenamebase = %s", filenamebase)

	# Continuous colormap
	Qt_filename = './' + filenamebase+'/Qdt_hist.csv'
	Qt_filename_std = './' + filenamebase+'/Qdt_std.csv'
	Qt = pd.read_csv(Qt_filename)
	Qt_std = pd.read_csv(Qt_filename_std)
	
	Tit_q = Tit + r', $\beta = $' + str(beta)
	savefig_name = f'./' + filenamebase + '/' + 'beta' + str(beta) + '_' + figname # To save plot
	
	qValues = np.arange(1, N+1)
	
	# setup the normalization and the colormap
	normalize = mcolors.Normalize(vmin = qValues.min(), vmax = qValues.max())
	colormap = cm.get_cmap('nipy_spectral')#, qValues.max()+1)
	
	Q_t = pl.figure()
	
	ax = pl.axes()
	ax = pl.axes(xscale = 'log', yscale = 'log')
	ax.grid()
	ax.set_xlim(Qt["t"].min(), Qt["t"].max())
	ax.set_ylim(1, N+2)
	pl.xlabel(r'$\ln(\delta t/f)$', fontsize = 23)
	pl.ylabel(r'$\ln\left(\langle Q(t) \rangle\right)$', fontsize = 23)
	pl.title(Tit_q, fontsize = 26, pad = 15)
	pl.xticks(fontsize = 16)
	pl.yticks(fontsize = 16)
	locs, labels = pl.yticks()
	pl.yticks(np.arange(1, N + 2, 2), np.arange(1, N + 2, 2))
	for q in qValues:
		#pl.scatter(Qt["t"], Qt["Q"+str(q)], color = colormap(normalize(q)))
		#ax.errorbar(Qt["t"], Qt["Q"+str(q)], yerr = Qt_std["Q"+str(q)], color = colormap(normalize(q)), linewidth = 1, ms = 1.5)
		pl.plot(Qt["t"], Qt["Q"+str(q)], color = colormap(normalize(q)), linewidth = 1, ms = 1.5)
		pl.fill_between(Qt["t"], Qt["Q"+str(q)] + Qt_std["Q"+str(q)], Qt["Q"+str(q)] - Qt_std["Q"+str(q)],#
			color = colormap(normalize(q)), alpha = 0.2, linewidth = 1)
	
	# setup the colorbar
	scalarmappaple = cm.ScalarMappable(norm = normalize, cmap = colormap)
	scalarmappaple.set_array(qValues)
	
	cbar = pl.colorbar(scalarmappaple, orientation = 'vertical', pad = 0.03, ticks = range(0, qValues.max()+1, 3))
	cbar.ax.set_xlabel(r'$q_i(t)$', fontsize = 23)
	
	pl.subplots_adjust(left = 0.075, bottom = 0.075, right = 1.05, top = 0.925)
	pl.savefig(savefig_name, dpi = 500)
	pl.close(Q_t)
	
	logging.info("Plots Q(t) for beta = %s created and saved in the root folder with name %s", beta, figname)















# Discrete colormap
'''Qt_filename = './' + filenamebase+'/Qdt_hist.csv'
Qt_filename_std = './' + filenamebase+'/Qdt_std.csv'
Qt = pd.read_csv(Qt_filename)
Qt_std = pd.read_csv(Qt_filename_std)

Tit_q = Tit + r', $\beta = $' + str(beta)
savefig_name = f'./' + filenamebase + '/' + 'beta' + str(beta) + '_' + figname # To save plot

qValues = np.arange(1, N+1)

hsv2rgb = lambda hue: mcolors.hsv_to_rgb([hue, 0.9, 0.7])
hues = np.linspace(0.75, 0.01, len(qValues))
colors = [hsv2rgb(hue) for hue in hues]

Q_t, ax = pl.subplots()

ax = pl.axes()
ax = pl.axes(xscale = 'log', yscale = 'log')
ax.grid()
ax.set_ylim(1, N+2)
ax.yaxis.set_minor_formatter(FormatStrFormatter("%.f"))
pl.xlabel(r'$\ln(\delta t/f)$', fontsize = 14)
pl.ylabel(r'$\ln\left(\langle Q(t) \rangle\right)$', fontsize = 14)
pl.title(Tit_q, fontsize = 16)
for q in qValues:
	#pl.scatter(Qt["t"], Qt["Q"+str(q)], color = colormap(normalize(q)))
	ax.errorbar(Qt["t"], Qt["Q"+str(q)], yerr = Qt_std["Q"+str(q)], color = colors[q-1], linewidth = 1, ms = 1.5)

# Fake a ScalarMappable so you can display a colormap
cmap, norm = mcolors.from_levels_and_colors(range(len(qValues) + 1), colors)
sm = cm.ScalarMappable(cmap = cmap, norm = norm)
sm.set_array([])

cbar = Q_t.colorbar(sm, orientation = 'vertical', pad = 0.02, ticks = range(1, N+1, 3))

cbar.set_label('$q_i(t)$', size = 16)
#cbar.set_ticks(np.arange(1, N+1, 4))

pl.savefig(savefig_name)
pl.close(Q_t)'''

























