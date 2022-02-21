import numpy as np
import pandas as pd
import glob
import argparse
import json
import logging
import pickle

def argument_parsing():
	""" Parse parameters filename for the analysis
	Returns: parser.parser_arg() : namespace
	Arguments from the command line """
	parser = argparse.ArgumentParser()
	parser.add_argument("filename", type=str, help="name of the JSON file with the simulations parameters")
	parser.add_argument("-s", "--steps", type=int, help="number of time steps to generate. Default = 100", default=100)
	parser.add_argument("--silent", help="do not print diagnose messages", action="store_true")
	parser.add_argument("--log",#
		help="sets the log level. Default=INFO. DEBUG will print detailed output of the simulation",#
		default="INFO")
	return parser.parse_args()

def get_params(filename):
	""" Read the parameters from the JSON file
	filename : str -> filename of the JSON file with simulation parameters data """
	with open(filename, "r") as read_file:
		params = json.load(read_file)
	return params

def build_filenamebase(args_dict):
	""" Builds the filenamebase.
	Paramemeters
		args : namespace with the program arguments
		Filenaming convention: 
		filename base is filename = Nxx_Niterxx_Nsimxx_Nsimstxx_betaxx_dtxx_lfxx_stfxx_ """
	filename = ''				# construct the filename with the simulation parameters:
	for key, val in sorted(args_dict.items()):
		if key != 'silent' and key != 'freediffusion' and key != 'Nf' and key != 'log':
			filename = filename + key + str(val)+'_'
		if args_dict['freediffusion'] == True:
			filename = filename + 'freediffusion_'
	return filename

def time_df(dt, N_iter, steps):
	""" Generates time dataframe (axis x) with steps uniformly on a log scale
    Parameters:
		dt: float			-> original time step
		N_iter: int			-> Total number of iterations (time steps) in the simulation
		steps: int			-> number of time steps to generate 
	Returns:
		t_df: pd.Dataframe	-> The time dataframe """
	tf = N_iter*dt
	t_df = pd.DataFrame([(t) for t in np.geomspace(dt, tf, steps)], columns = ['t'])
	return t_df

def empty_Q_hist(N_iter, dt, Qmax, steps = 100):
	""" Builds an empty dataframe for the histograms of charges with time steps uniformly distributed on a log scale
	Parameters:
		N_iter: int		-> Total number of time steps in the simulation
		dt: float		-> time step of the simulation
		Qmax: int		-> total number posible charges per particle
		steps : int		-> Number of time steps to generate
		Default steps = 100
	Returns:
		hists_out: pd.DataFrame -> Empty data frame for charge histograms"""
	hists_out = time_df(dt, N_iter, steps)
	for Q in range(1, Qmax + 1):
		col_name = 'Q' + str(Q)
		hists_out[col_name] = 0.0
	return hists_out

def empty_Q_pd(N_iter, Nsims, dt, Qmax, steps = 100):
	""" Builds an empty dataframe for the histograms of charges with time steps uniformly distributed on a log scale
	Parameters:
		N_iter: int		-> Total number of time steps in the simulation
		dt: float		-> time step of the simulation
		Qmax: int		-> total number posible charges per particle
		steps : int		-> Number of time steps to generate
		Default steps = 100
	Returns:
		hists_out: pd.DataFrame -> Empty data frame for charge histograms"""
	hists_out = time_df(dt, N_iter, steps)
	for Q in range(1, Qmax + 1):
		for sim in range(1, Nsims + 1):
			col_name = 'Sim' + str(sim) + 'Q' + str(Q)
			hists_out[col_name] = 0.0
	return hists_out

def counters(filename, N, Nsims):
	counts = []
	for sim in range(Nsims):
		N_min = N
		counter = 0
		sim_col = 'N_' + str(sim+1)
		for i in filename.index:
			if filename.at[i, sim_col] < N_min:
				N_min = filename.at[i, sim_col]
				counter = i
		counts.append(counter)
	return counts

def decompress_data(compressed_df, N, N_iter, Nsims, dt, steps = 100):
	""" Decompresses a data frame that has only the time when the particles changes to one with several times uniformly
		spaced on a log scale.
	Parameters:
		compressed_df: pd.DataFrame -> compressed dataframe with only the changes in N and the times
		dt: float		-> original time step
		N_iter: int		-> Total number of iterations (time steps) in the simulation
		steps : int		-> number of time steps to generate
	Returns:
		pd.DataFrame	-> decompressed data with all the requested times repeating the number of particles when it does not change"""
	counts = counters(compressed_df, N, Nsims)
	decomp_df = time_df(dt, N_iter, steps)
	for sim in range(Nsims):
		tp = counts[sim]+1 # Total sample points for each sim in data DT
		N_col = 'N_' + str(sim+1)
		t_col = 't_' + str(sim+1)
		for i in range(tp):
			t_now = float(compressed_df.at[i, t_col])
			if(i == tp-1):
				tf = N_iter*dt
			else:
				tf = float(compressed_df.at[i+1, t_col])
			N_now = float(compressed_df.at[i, N_col])
			decomp_df.loc[(t_now <= decomp_df['t']) & (decomp_df['t'] <= tf), N_col] = N_now
	return decomp_df

def decompress_data_Q(compressed_hists, N, N_iter, Nsims, dt, steps = 100):
	""" Decompresses an array that has only the time when the particle charge histogram changes to one with many times uniformly distributed on a log scale.
	Parameters:
		N: int -> Initial number of particles
		compressed_hists: array -> Compressed array of histograms with only the changes in the charges and the times
		dt: float -> original time step
		N_iter: int -> Total number of iterations (time steps) in the simulation
		steps: int -> Number of time steps to generate. 
		Default steps = 100
	Returns:
		hists_out: pd.Dataframe of histograms -> decompressed data with all the times repeating the histograms do not change"""
	hists_out = empty_Q_pd(N_iter, Nsims, dt, N, steps)
	for sim in range(Nsims):
		for i, line in enumerate(compressed_hists[sim]):  #range(len(hists[Sim])):
			s, t_now, (hist, bins) = compressed_hists[sim][i] #[i][1]
			if i < len(compressed_hists[sim]) - 1:
				tf = compressed_hists[sim][i + 1][1]
			else:
				tf = N_iter*dt
			for idx, Q in enumerate(hist):
				col_name = 'Sim' + str(sim + 1) + 'Q' + str(idx + 1)
				hists_out.loc[(t_now <= hists_out['t']) & (hists_out['t'] <= tf), col_name] = Q
	return hists_out


























        






def build_finalfilenames(filenamebase, Nsims, NsimsQ, steps):
    """ Filenames to save the data analysis

    Parameters:
    ===========
    filenamebase : str
        filename base obtained from build_filenames(params_dict)
    Nsims : int
        Total number of simulations (number CVS files for N evolution found)
    NsimsQ : int
        Total number of simulations (number pickle files for charge evolution found)

    Returns
    filenamefinal : str
        filename (CSV) to save data for number of particles N evolution
    filenamefinalQ :str
        filename (CSV) to save data for number of charge Q evolution
    """

    filenamefinal = 'resumen_N_'+filenamebase+'Nsimstot'+str(Nsims)+'_steps'+str(steps)+'.csv'
    filenamefinalQ = 'resumen_Q_'+filenamebase+'Nsimstot'+str(NsimsQ)+'_steps'+str(steps)+'.csv'
    return filenamefinal, filenamefinalQ





def analyse_Q(Qmax, N_iter, dt, steps, filenames):
    """ 
    Builds the dataframe with the evolution of the number of particles:
    t log_t N_avg std_N log_N

    Parameters:
    ===========
    Qmax : int
        Maximum charge for a particle
    N_iter : int
        Total number of time steps of each simulation
    dt : float
        time step of each simulation
    steps : int
        Number of steps to generate
    filenames : list 
        list of all the pickle filenames to analyse with the evolution of the histograms 
        of charge for each simulation

    Returns:
    ========
    Q_df : pd.DataFrame
        Dataframe with the evolution of histogram charges
        Columns:
        t, Q1, Q2, ..., Q100 
    """
    # initial charge of the particles
    Q_df=empty_Q_pd(N_iter, dt, Qmax, steps)
    count=1
    Nsims=len(filenames)
    for filename in filenames:
        logging.info("Analysing pickle file # %s : %s", count, filename)
        # load the histogram file
        with open(filename, 'rb') as histsfile:
            hists=pickle.load(histsfile)
        tmp_df=decompress_data_Q(N,hists,dt,N_iter,steps)
        # do the average ... 
        for col in Q_df:
            if(col != 't'):
                Q_df[col]=Q_df[col]+tmp_df[col]    
        count=count+1
    for col in Q_df:
        if(col != 't'):
            Q_df[col]=Q_df[col]/Nsims 
    return Q_df
