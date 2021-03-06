import numpy as np
import pandas as pd
import os
import glob
import argparse
import json
import logging
import pickle
import matplotlib.pyplot as pl
import matplotlib.ticker as mticker
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
pl.style.use('classic')

import Functions as Fun

args = Fun.argument_parsing()
loglevel = args.log
if not args.silent :
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(format='%(asctime)s : %(message)s', level=numeric_level)

params = Fun.get_params(args.filename)
logging.info("Analysis program called with arguments:\n%s", json.dumps(vars(args), indent=4, sort_keys=True))
logging.info("Starting analysis using the parameters:\n %s", json.dumps(params, indent=4, sort_keys=True))

filenamebase = Fun.build_filenamebase(params)
logging.debug("filenamebase = %s", filenamebase)

N = int(params["N"])
N_iter = int(params["Niter"])
Nsims = int(params["Nsim"])
beta = float(params["beta"])
dt = float(params["dt"])
lf = float(params["lf"])

Qmax = N
steps = args.steps

logging.info("Starting read data and statistical analysis for N(t)")
Nt_filename = './' + filenamebase+'/Nsim_' + str(Nsims) + '_N_t.csv'
DT = pd.read_csv(Nt_filename)

decomp_df = Fun.decompress_data(DT, N, N_iter, Nsims, dt, steps)

N_df = pd.DataFrame(columns = ['t', 'N_avg', 'std_N', 'log_t', 'log_N'])
N_df['t'] = decomp_df['t']
N_df['N_avg'] = decomp_df.iloc[:, 1:].mean(axis = 1)
N_df['std_N'] = decomp_df.iloc[:, 1:].std(axis = 1)
N_df['log_t'] = np.log(decomp_df['t'])
N_df['log_N'] = np.log(N_df['N_avg'])
logging.info("Finishing read data and statistical analysis for N(t)")

Ndf_filename = './' + filenamebase+'/Ndt.csv'
N_df.to_csv(Ndf_filename, index = False)
logging.info("Wrote particle number evolution data analysis to %s", Ndf_filename)


logging.info("Starting read data and statistical analysis for Q(t) histograms")
hist_filename = './' + filenamebase+'/sim_' + str(Nsims) + '_Q_t.pickle'
with open(hist_filename, 'rb') as histsfile:
	hists=pickle.load(histsfile)

tmp_df = Fun.decompress_data_Q(hists, N, N_iter, Nsims, dt, steps)
Q_hist_ave = Fun.empty_Q_hist(N_iter, dt, Qmax, steps)
Q_hist_std = Fun.empty_Q_hist(N_iter, dt, Qmax, steps)

for q in range(1, Qmax + 1):
	col_name = 'Q' + str(q)
	Q_hist_ave[col_name] = tmp_df.loc[:, 'Sim1Q' + str(q):'Sim'+str(Nsims)+'Q' + str(q)].mean(axis = 1)
	Q_hist_std[col_name] = tmp_df.loc[:, 'Sim1Q' + str(q):'Sim'+str(Nsims)+'Q' + str(q)].std(axis = 1)	
logging.info("Finishing read data and statistical analysis for Q(t) histograms")

Q_hist_filename = './' + filenamebase+'/Qdt_hist.csv'
Q_hist_ave.to_csv(Q_hist_filename, index = False)

Q_hist_std_filename = './' + filenamebase+'/Qdt_std.csv'
Q_hist_std.to_csv(Q_hist_std_filename, index = False)
logging.info("Wrote charge histogram evolution data analysis to %s", Q_hist_filename)

'''Tit = r'$N_0 = $' + str(N) + r', $l_f = $' + str(lf)
lab = r'$\beta = $' + str(beta)
figname = 'Ninit' + str(N) + '_lf' + str(lf) + '_dt' + str(dt) + '.pdf'
savefig_name = f'./' + filenamebase + '/' + figname # To save plot

N_avg_t = pl.figure()
ax = pl.axes()
ax = pl.axes(xscale = 'log', yscale = 'log')
ax.grid()
ax.set_ylim(40, N+2)
#ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_formatter(FormatStrFormatter("%.f"))
pl.xlabel(r'Time $(\delta t/f)$', fontsize = 14)
pl.ylabel(r'$\langle N(t) \rangle$', fontsize = 14)
pl.title(Tit, fontsize = 16)
ax.errorbar(N_df['t'], N_df['N_avg'], yerr = N_df['std_N'], label = lab, linewidth = 0.75, ms = 1.5)
pl.legend(loc = 'lower left', fontsize = 10)
N_avg_t.savefig(savefig_name, dpi = 500)
pl.close(N_avg_t)'''











































