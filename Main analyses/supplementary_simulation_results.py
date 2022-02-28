## simulate differences in HRF shape

import numpy as np
import matplotlib.pyplot as plt
from Simulate_data import Simulations
from deconvolution_tools import data_deconvolution_canon
from gsbs import GSBS
from time_correlation_simple import time_correlation_simple

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

basedir = '/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/'
savedir_final = basedir + 'results_fast_event_complete/'

#default settings
nvox = 50
ntime = 200
nstates = 15
nsub = 1
group_std = 0
sub_std = 1
TR = 2.47
reps = 100
length_std = 1
maxK = 100
peak_delay = 6
peak_disp = 1
extime = 2
maxK=100
nstates_max = 100

nstates = np.full(3, 0).astype(int)
tdists = np.full([3, maxK + 2], 0).astype(float)
bounds = np.full([3, ntime], 0).astype(int)

#run all the simulations
sim=Simulations(nvox=nvox, ntime=ntime, nstates=50, nsub=nsub, group_std=group_std, TR=TR, length_std=length_std,peak_delay=peak_delay,peak_disp=1, extime=extime, sub_std=sub_std, maxK=maxK)
bounds50,data50, subbounds50  = sim.generate_simulated_data_HRF(rep=1)

sim=Simulations(nvox=nvox, ntime=ntime, nstates=15, nsub=nsub, group_std=group_std, TR=TR, length_std=length_std,peak_delay=peak_delay,peak_disp=1, extime=extime, sub_std=sub_std, maxK=maxK)
bounds15,data15, subbounds15  = sim.generate_simulated_data_HRF(rep=1)

sim=Simulations(nvox=nvox, ntime=ntime, nstates=30, nsub=nsub, group_std=group_std, TR=TR, length_std=length_std,peak_delay=peak_delay,peak_disp=1, extime=extime, sub_std=sub_std, maxK=maxK)
bounds30,data30, subbounds30  = sim.generate_simulated_data_HRF(rep=1)

sim=Simulations(nvox=nvox, ntime=ntime, nstates=50, nsub=nsub, group_std=group_std, TR=TR, length_std=length_std,peak_delay=5.5,peak_disp=1, extime=extime, sub_std=sub_std, maxK=maxK)
bounds50_55,data50_55, subbounds50_55 = sim.generate_simulated_data_HRF(rep=1)

sim=Simulations(nvox=nvox, ntime=ntime, nstates=50, nsub=nsub, group_std=group_std, TR=TR, length_std=length_std,peak_delay=4.6,peak_disp=1, extime=extime, sub_std=sub_std, maxK=maxK)
bounds50_46,data50_46, subbounds50_46  = sim.generate_simulated_data_HRF(rep=1)

## show how the number of states are underestimated in this case
plot_states=15
fig, axs = plt.subplots(2,3)
StateSegment = GSBS(x=np.squeeze(data15), kmax=maxK, statewise_detection=True)
StateSegment.fit()
tdists[0, :] = StateSegment.tdists
nstates[0] = StateSegment.nstates
bounds[0, :] = StateSegment.bounds
time_correlation_simple(axs[0,0], np.corrcoef(StateSegment.x),StateSegment.bounds)
axs[1,0].plot(StateSegment.tdists)
axs[1,0].vlines(x=plot_states, ymin=0, ymax=StateSegment.tdists.max())
axs[1,0].set_title('15 states')

plot_states=30
StateSegment = GSBS(x=np.squeeze(data30), kmax=maxK, statewise_detection=True)
StateSegment.fit()
tdists[0, :] = StateSegment.tdists
nstates[0] = StateSegment.nstates
bounds[0, :] = StateSegment.bounds
time_correlation_simple(axs[0,1], np.corrcoef(StateSegment.x),StateSegment.bounds)
axs[1,1].plot(StateSegment.tdists)
axs[1,1].vlines(x=plot_states, ymin=0, ymax=StateSegment.tdists.max())
axs[1,1].set_title('30 states')

plot_states=50
StateSegment = GSBS(x=np.squeeze(data50), kmax=maxK, statewise_detection=True)
StateSegment.fit()
tdists[0, :] = StateSegment.tdists
nstates[0] = StateSegment.nstates
bounds[0, :] = StateSegment.bounds
time_correlation_simple(axs[0,2], np.corrcoef(StateSegment.x),StateSegment.bounds)
axs[1,2].plot(StateSegment.tdists)
axs[1,2].vlines(x=plot_states, ymin=0, ymax=StateSegment.tdists.max())
axs[1,2].set_title('50 states')
plt.tight_layout()
plt.savefig(savedir_final + 'anticorrelated_states_K.pdf')

## show that estimating the HRF incorrectly does not have a big impact
plot_states=50
fig, axs = plt.subplots(2,3)
data = data_deconvolution_canon(np.squeeze(data50), TR)
StateSegment = GSBS(x=data, kmax=maxK,  statewise_detection=True)
StateSegment.fit()
tdists[0, :] = StateSegment.tdists
nstates[0] = StateSegment.nstates
bounds[0, :] = StateSegment.bounds
time_correlation_simple(axs[0,0], np.corrcoef(StateSegment.x),StateSegment.bounds)
axs[1,0].plot(StateSegment.tdists)
axs[1,0].vlines(x=plot_states, ymin=0, ymax=StateSegment.tdists.max())
axs[1,0].set_title('HRF 5')

data = data_deconvolution_canon(np.squeeze(data50_46), TR)
StateSegment = GSBS(x=data, kmax=maxK,  statewise_detection=True)
StateSegment.fit()
tdists[1, :] = StateSegment.tdists
nstates[1] = StateSegment.nstates
bounds[1, :] = StateSegment.bounds
time_correlation_simple(axs[0,1], np.corrcoef(StateSegment.x), StateSegment.bounds)
axs[1,1].plot(StateSegment.tdists)
axs[1,1].vlines(x=plot_states, ymin=0, ymax=StateSegment.tdists.max())
axs[1,1].set_title('HRF 4.6')

data = data_deconvolution_canon(np.squeeze(data50_55), TR)
StateSegment = GSBS(x=data, kmax=maxK,  statewise_detection=True)
StateSegment.fit()
tdists[2, :] = StateSegment.tdists
nstates[2] = StateSegment.nstates
bounds[2, :] = StateSegment.bounds
time_correlation_simple(axs[0,2], np.corrcoef(StateSegment.x), StateSegment.bounds)
axs[1,2].plot(StateSegment.tdists)
axs[1,2].vlines(x=plot_states, ymin=0, ymax=StateSegment.tdists.max())
axs[1,2].set_title('HRF 5.5')
plt.tight_layout()

plt.savefig(savedir_final + 'HRF_differences.pdf')
