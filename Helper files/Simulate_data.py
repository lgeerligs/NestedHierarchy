import numpy as np
from hrf_estimation import hrf
savedir = '/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/simulations/'

class Simulations:
    def __init__(self, nvox, ntime, nstates, nsub, sub_std, group_std, TR, length_std,maxK, peak_delay:float=6, peak_disp:float=1, extime:int=2):
        self.nvox=nvox
        self.ntime=ntime
        self.nstates=nstates
        self.nsub=nsub
        self.group_std=group_std
        self.sub_std=sub_std
        self.TR=TR
        self.length_std=length_std
        self.peak_delay=peak_delay
        self.peak_disp=peak_disp
        self.extime = extime
        self.maxK=maxK


    @staticmethod
    def sample(n: int, p: float, r: float):
        # n = number of timepoints
        # p = number of states
        # r = variability of state lengths
        x_0 = np.linspace(0, n, p + 1).astype(int)[1: -1]
        x_r = []
        q = r * (n / p) - 0.5

        for i in range(p - 1):
            while True:
                x = x_0[i] + np.random.randint(-q, q + 1)
                if (i > 0 and x_r[i - 1] >= x) or (i == 0 and x < 1) or (x > n-(p-i)-1):
                    continue
                else:
                    break

            x_r.append(x)

        x_r = np.array(x_r)
        # x_d = np.concatenate(([x_r[0]], x_r[1:] - x_r[:-1], [n - x_r[-1]]))
        bounds = np.zeros(n).astype(int)
        bounds[x_r]=1
        states = Simulations.deltas_states((bounds))
        print(max(states))
        # x_d = np.concatenate(([x_r[0]], x_r[1:] - x_r[:-1], [n - x_r[-1]]))

        return bounds, states #x_d, min(x_d), max(x_d)

    def generate_simulated_data_HRF(self, nstates=None, group_std=None, sub_std=None, length_std=None, peak_delay=None, peak_disp=None, extime=None, TR=None, TRfactor=1, nsub=None, rep=500):

        nstates = nstates or self.nstates
        group_std = group_std or self.group_std
        sub_std = sub_std or self.sub_std
        length_std = length_std or self.length_std
        peak_delay = peak_delay or self.peak_delay
        peak_disp = peak_disp or self.peak_disp
        extime = extime or np.int(self.extime/TRfactor)
        TR = TR or self.TR*TRfactor
        nsub = nsub or self.nsub

        np.random.seed(rep)
        state_means = np.random.randn(self.nvox, nstates)

        bounds,state_labels = Simulations.sample(self.ntime, nstates, length_std)

        evData = np.zeros([self.ntime+extime, self.nvox])
        for t in range(0,self.ntime):
            if (state_labels[t] % 2) == 0:
                evData[t,:] = state_means[:, 0]
            else:
                evData[t, :] = state_means[:, 0] * -1

        #extend the final state to make sure it is present in the final signal
        for te in range(self.ntime, self.ntime+extime):
            evData[te,:]=evData[-1,:]

        spmhrf = hrf.spm_hrf_compat(np.arange(0,30,TR), peak_delay=peak_delay, peak_disp=peak_disp)
        BOLDevData = self.convolve_with_hrf(evData, spmhrf, extime)

        groupSig = group_std * np.random.randn(self.ntime+extime, self.nvox)
        BOLDgroupSig = self.convolve_with_hrf(groupSig, spmhrf, extime)

        subData = np.zeros([nsub,  self.ntime, self.nvox])
        subbounds = np.zeros([nsub, self.ntime])
        for s in range(0,nsub):
            subSig = sub_std * np.random.randn(self.ntime + extime, self.nvox)
            subSig=sub_std * subSig
            BOLDsubSig = self.convolve_with_hrf(subSig, spmhrf, extime)
            subData[s,:,:] = BOLDsubSig + BOLDgroupSig + BOLDevData


        return bounds, subData, subbounds

    def convolve_with_hrf(self, signal, hrf, extime):
        BOLDsignal = np.zeros([self.ntime + extime + len(hrf) - 1, self.nvox])
        for n in range(0, self.nvox):
            BOLDsignal[:,n] = np.convolve(signal[:,n], hrf)
        BOLDsignal = BOLDsignal[extime:extime + self.ntime,:]
        return BOLDsignal

    @staticmethod
    def deltas_states(deltas: np.ndarray) -> np.ndarray:
        deltas.astype(int)
        states = np.zeros(deltas.shape[0], int)
        for i, delta in enumerate(deltas[1:]):
            states[i + 1] = states[i] + 1 if delta else states[i]

        return states

