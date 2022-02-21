# Compile directly from root source in terminal.
# To compile the cython helpers modules use in terminal:
#	make all
# Then:
#	python3 Dyson_Fusion_Gas.py 20 1000 1 1 0.5 0.00001 0.5 0
#	sintax:
#		python3 Sim.py N Niter Nsim Nsimst beta dt lf stf
'''
	N				->	Initial number of particles
	Niter		->	Iteration number per simulation
	Nsim		->	Number of simulations to perform
	Nsimst	->	Number label for first simulation
	beta		->	Inverse temperature
	dt			->	Virtual time step
	lf			->	fusion length
	stf			->	Iteration number to start fusion process (To include a time for thermal equilibrium)
'''
# These files need to be in the same directory

import numpy as np
import pandas as pd
import os
import argparse
import pickle
import matplotlib.pyplot as plt
from operator import itemgetter, attrgetter
import realfusionhelpers as rfh

import json
import logging
import time
import Functions as fun

''' ***********************************************Main program****************************************************************** '''
inicio = time.time()
args = fun.argument_parsing()

	# Set up logging info assuming loglevel is bound to the string value obtained from the command line argument.
	# Convert to upper case to allow the user to specify --log=DEBUG or --log=debug
loglevel = args.log
if not args.silent:
	numeric_level = getattr(logging, loglevel.upper(), None)
	if not isinstance(numeric_level, int):
		raise ValueError('Invalid log level: %s' % loglevel)
	logging.basicConfig(format='%(asctime)s : %(message)s', level=numeric_level)

filenamebase = fun.build_filenamebase(args)
	# Directory to save all data:
os.makedirs(filenamebase, exist_ok = True)
fun.saveargs(args, './'+filenamebase+'/')

	# Run the simulations:
logging.info("Starting %s simulations with parameters:\n %s", args.Nsim, json.dumps(vars(args), indent=4, sort_keys=True))

	# Time-evolution of frequency charges data:
Q_evol = []
	# Time evolution of N for each simulation. Allows make the statistics:
N_evol = pd.DataFrame()

for sim in range(args.Nsimst, args.Nsimst + args.Nsim):
		# Makes one simulation only:
	data_Ntime, data_Qtime = fun.one_simulation_no_file(args.N, args.Niter, args.dt, args.beta, args.lf,#
			args.stf, sim, args.freediffusion, args.Nf)

		# Time evolution of number particles in the system for each value of sim:
	N_evol = pd.concat([N_evol, data_Ntime], axis = 1)
		# Time evolution of the system's total charge for each value of sim:
	Q_evol.append(data_Qtime)

	# Save the data on a file:
filename = './'+filenamebase+'/Nsim_'+str(sim)+'_N_t.csv'
N_evol.to_csv(filename, index = False)
filenamehists = './'+filenamebase+'/Nsim_'+str(sim)+'_Q_t.pickle'
with open(filenamehists, 'wb') as histsfile:
	pickle.dump(Q_evol, histsfile)
logging.debug("%s", Q_evol)
logging.info("Wrote particle number results to %s", filename)
logging.info("Wrote charges histogram results to %s", filenamehists)

logging.info("Finished %s simulations", args.Nsim)

fin = time.time()
print('Total time for make all simulations ' + str(fin-inicio))
''' ***************************************************************************************************************************** '''


'''
for sim in range(args.Nsimst, args.Nsimst + args.Nsim):

	#DT_Sim = pd.DataFrame(columns = ['time', 'label', 'charge', 'pos', 'step', 'newpos'])	# Dataframe for particles evolution
	#DT_Sim_Name = './'+ filenamebase + '/Sim'+str(sim)+'_Data.csv'												# Name of file to save data
	
	#DT_Sim_File = open(DT_Sim_Name, 'w')	# Open the file to write all evolution of each particle and each sim
	data_Ntime, data_Qtime = fun.one_simulation_no_file(#
			args.N, args.Niter, args.dt, args.beta, args.lf,#
			args.stf, sim, args.freediffusion, args.Nf)		# Makes one simulation only
	#DT_Sim_File.close()						# Close the file with all evolution of each particle and each sim data
	
	N_evol = pd.concat([N_evol, data_Ntime], axis = 1)		# Time evolution of number particles in the system for each value of sim
	Q_evol.append(data_Qtime)								# Time evolution of the system's total charge for each value of sim
'''



























