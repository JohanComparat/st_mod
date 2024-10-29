
#import matplotlib
#matplotlib.use('Agg')
#matplotlib.rcParams.update({'font.size': 14})
#import matplotlib.pyplot as plt
import time
t0 = time.time()
import os, sys
import glob
import numpy as np

import ultranest
from ultranest.plot import cornerplot
from ultranest.plot import PredictionBand
from scipy.special import erf

dir_2_OUT = 'fitDATAlxms'
params0 = np.array([20, 2, 11, 0.1, 4, 2.])
params1 = np.array([41.1, 1.8, 11.25, 0.1, 4.1])

def fun_Ms_Mh_LX(log10Ms, params):
	log10Mh_0, alpha1, log10M0, delta, alpha2, alpha3 = params
	log10Mh_part1 = log10Mh_0 + alpha1 * (log10Ms - log10M0)
	log10Mh_part2 = delta * ( alpha2 - alpha1) * np.log10( 0.5 * (1. + 10**(log10Ms - log10M0) )**(1./delta) )
	log10LX = alpha3 * (log10Mh_part1 + log10Mh_part2)
	return log10LX

def fun_Ms_LX(log10Ms, params):
	log10Mh_0, alpha1, log10M0, delta, alpha2 = params
	log10Mh_part1 = log10Mh_0 + alpha1 * (log10Ms - log10M0)
	log10Mh_part2 = delta * ( alpha2 - alpha1) * np.log10( 0.5 * (1. + 10**(log10Ms - log10M0) )**(1./delta) )
	log10LX = log10Mh_part1 + log10Mh_part2
	return log10LX

x_data = np.arange(9.5, 12, 0.01)
y_data = fun_Ms_LX(x_data, params1)
y_data_err = 0.01

x_data, y_data, x_data_err, y_data_err = np.loadtxt('LxMs.ascii', unpack = True)

def fun_Ms_LX_mod(params):
	log10Mh_0, alpha1, log10M0, delta, alpha2 = params
	log10Mh_part1 = log10Mh_0 + alpha1 * (x_data - log10M0)
	log10Mh_part2 = delta * ( alpha2 - alpha1) * np.log10( 0.5 * (1. + 10**(x_data - log10M0) )**(1./delta) )
	log10LX = log10Mh_part1 + log10Mh_part2
	return log10LX

param_names = [r'$\log_{10}(M_{h0})$', r'$\alpha_{1}$', r'$\log_{10}(M_0)$', r'$\delta$', r'$\alpha_{2}$', r'$\alpha_{3}$']

def my_prior_transform(cube):
	"""
	returns a cube of parameter boundaries
	log10Mh_0, alpha1, log10M0, delta, alpha2, alpha3
	"""
	params = cube.copy()
	# p0 log10Mh_0
	lo =  18.0
	hi = 22.0
	params[0] = cube[0] * (hi - lo) + lo
	# p1 alpha1
	lo = 1
	hi = 3
	params[1] = cube[1] * (hi - lo) + lo
	# p2 log10M0
	lo = 10
	hi = 12
	params[2] = cube[2] * (hi - lo) + lo
	# p1 delta
	lo = 0.
	hi = 0.2
	params[3] = cube[3] * (hi - lo) + lo
	# p1 alpha2
	lo = 3
	hi = 5
	params[4] = cube[4] * (hi - lo) + lo
	# p1 alpha3
	lo = 1
	hi = 3
	params[5] = cube[5] * (hi - lo) + lo
	return params

param_names = [r'$\log_{10}(L_{0})$', r'$\alpha_{1}$', r'$\log_{10}(M_0)$', r'$\delta$', r'$\alpha_{2}$']#, r'$\alpha_{3}$']

def my_prior_transform(cube):
	"""
	returns a cube of parameter boundaries
	log10Mh_0, alpha1, log10M0, delta, alpha2, alpha3
	"""
	params = cube.copy()
	# p0 log10Mh_0
	lo =  40.0
	hi = 42.0
	params[0] = cube[0] * (hi - lo) + lo
	# p1 alpha1
	lo = 1
	hi = 3
	params[1] = cube[1] * (hi - lo) + lo
	# p2 log10M0
	lo = 10
	hi = 12
	params[2] = cube[2] * (hi - lo) + lo
	# p1 delta
	lo = 0.
	hi = 0.2
	params[3] = cube[3] * (hi - lo) + lo
	# p1 alpha2
	lo = 3
	hi = 5
	params[4] = cube[4] * (hi - lo) + lo
	return params

def my_likelihood(params):
	"""
	returns the likelihood value
	"""
	y_model = fun_Ms_LX_mod(params)
	like = np.sum( -0.5 * ( (y_model - y_data)/y_data_err )**2 )
	return like

print('initializes sampler', time.time() - t0)
sampler = ultranest.ReactiveNestedSampler(
	param_names,
	loglike=my_likelihood,
	transform=my_prior_transform,
	log_dir = dir_2_OUT,
	resume='overwrite' #or 'resume' or 'overwrite' or 'subfolder',
	)

print('run fitting', time.time() - t0)
result = sampler.run(
	min_num_live_points = 400  ,            # 400,
	dlogz = 0.5,                          # dlogz=0.5, # desired accuracy on logz
	# min_ess = 400.,                       # min_ess=400, # number of effective samples
	# update_interval_iter_fraction = 0.4, # how often to update region
	# max_num_improvement_loops = 2,       # how many times to go back and improve
	# cluster_num_live_points=40
	)

sampler.print_results()
print('finished fitting the parameters', time.time()-t0)
print('starts tabulating outputs')
print( '-' * 100 )
#q = 0.341
#q2 = 0.45
#q3 = 0.49


