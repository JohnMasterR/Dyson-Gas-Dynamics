# python3 Main_Dyson_Gas.py

import numpy as np
import pandas as pd
import os
import argparse
import pickle
import matplotlib.pyplot as plt
from operator import itemgetter, attrgetter

import logging
import time

inicio = time.time()


''' Simulation parameters: '''
N = 50
Niter = 200
Nsim = 10
Nsimst = 1
#dt = 1e-05
lf = 0.5
stf = 0
steps = 100

beta_min = 0.2
beta_max = 0.2
delta_beta = 0.2

lf_max = 0.5
delta_lf = 0.5

betas = np.arange(beta_min, beta_max+delta_beta, delta_beta)

base = 'N'+str(N)+'_Niter'+str(Niter)+'_Nsim'+str(Nsim)+'_Nsimst'+str(Nsimst)

while lf <= lf_max:
	for beta in betas:
		dt = 0.02*beta*(0.1*np.pi*lf/N)*(0.1*np.pi*lf/N)
		print('------ Working with ------ lf = '+str(lf)+' beta = '+str(beta))
		params = str(N)+' '+str(Niter)+' '+str(Nsim)+' '+str(Nsimst)+' '+str(beta)+' '+str(dt)+' '+str(lf)+' '+str(stf)
		file = base+'_beta'+str(beta)+'_dt'+str(dt)+'_lf'+str(lf)+'_stf'+str(stf)+'_/args.json'
		
			# Makes all data for each simulation:
		os.system('python3 ./Data_Script_Generator/Dyson_Fusion_Gas.py '+params)
			# Males statistical analysis for all sims:
		os.system('python3 ./Stat_Analysis/N_Q_Stat_Analysis.py '+file+' --steps '+str(steps))
		
	params = str(N)+' '+str(Niter)+' '+str(Nsim)+' '+str(Nsimst)+' '+str(dt)+' '+str(lf)+' '+str(stf)+' '+str(betas[0])+' '+str(betas[len(betas)-1])+' '+str(delta_beta)
	os.system('python3 ./Plots/Plots_Nt_Qt.py '+params)
	
	lf = lf + delta_lf


fin = time.time()
print('Total time for make all process ' + str(fin-inicio))

















