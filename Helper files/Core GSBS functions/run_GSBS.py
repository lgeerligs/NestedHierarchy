from gsbs import GSBS
import os
import numpy as np
import h5py
import pickle
from scipy.io import loadmat
from deconvolution_tools import data_deconvolution_canon
from gsbs_legacy import GSBS as GSBS_legacy

# run event detection on fMRI data
def run_state_detection_fixgroup(r:int, kfold_data:int, maxk:int, savename:str, savedir:str, basedir:str, statewise:bool=True, overwrite:bool=False, deconv:str='none', TR:float=0, data_denoise:bool=False, legacy_GSBS:bool=False):

    filename = savedir + 'GSevents_sl' + str(r) + '_kfold_data' + str(kfold_data)  + '_maxk' + str(maxk) + '_' + savename + '.pickle'

    if overwrite is True or not os.path.isfile(filename):
        #load data
        data, hr = loaddata(basedir, r, kfold_data, 0, data_denoise, deconv, TR)
        Ns = np.shape(data)[0]
        Nvox = np.shape(data)[1]

        if Nvox > 14:
            print(r)
            data = np.zeros((kfold_data, Ns, Nvox))
            hrf_shape = np.zeros((kfold_data,3))

            for k in range(0, kfold_data):
                   data[k, :, :], hrf_shape[k,:] = loaddata(basedir, r, kfold_data, k, data_denoise, deconv, TR)

            nstates = np.full(kfold_data, 0).astype(int)
            tdists = np.full([kfold_data, maxk + 2], 0).astype(float)
            bounds = np.full([kfold_data, Ns], 0).astype(int)
            strengths = np.full([kfold_data, Ns], 0).astype(float)
            meansig = np.full([kfold_data, Ns], 0).astype(float)
            bounds19 = np.full([kfold_data, Ns], 0).astype(int)
            strengths19 = np.full([kfold_data, Ns], 0).astype(float)

            for kf in np.arange(0,kfold_data):
                    if legacy_GSBS:
                        StateSegment = GSBS_legacy(x=np.squeeze(data[kf, :, :]), kmax=maxk)
                    else:
                        StateSegment = GSBS(x=np.squeeze(data[kf,:,:]), kmax=maxk, statewise_detection=statewise)

                    StateSegment.fit()

                    if legacy_GSBS:
                        tdists[kf, 0:-1] = StateSegment.tdists
                    else:
                        tdists[kf, :] = StateSegment.tdists
                    nstates[kf] = StateSegment.nstates
                    bounds[kf, :] = StateSegment.bounds
                    strengths[kf, :] = StateSegment.strengths
                    if statewise == False:
                        bounds19[kf, :] = StateSegment.get_bounds(k=20)
                        strengths19[kf, :] = StateSegment.get_strengths(k=20)
                    else:
                        bounds19[kf, :] = StateSegment.get_bounds(k=19)
                        strengths19[kf, :] = StateSegment.get_strengths(k=19)
                    meansig[kf, :] = np.mean(data[kf,:,:], axis=1)

            results = {'nstates': nstates, 'tdists': tdists, 'bounds': bounds, 'strengths': strengths,
                       'meansig': meansig, 'bounds19': bounds19, 'strengths19': strengths19, 'hrf_shape':hrf_shape}

            with open(filename, 'wb') as output:
                pickle.dump(results, output, pickle.HIGHEST_PROTOCOL)

            return results


def loaddata(basedir, r, kfold_data, curkfold, data_denoise, deconv, TR):

    if deconv == 'hrf_est':
        d = loadmat(basedir + 'data/ROI_' + str(r) + str(kfold_data) + 'group_hyp_step2_radius_3_step_2_HOA_cortex_group' + str(curkfold + 1) + 'hrf_denoised_data.mat')
        if data_denoise:
            data = d['data_denoised']
        else:
            data = d['data_notdenoised']

        hrf_shape = d['hrf_shape_all'].mean(axis=2).mean(axis=1)

    else:
        hrf_shape = np.zeros(3)
        if data_denoise == False and deconv == 'none':
            f = h5py.File(basedir + 'data/ROI_' + str(r) + str(kfold_data) + 'group_hyp_step2_radius_3_step_2_HOA_cortex_group' + str(curkfold + 1) +  '.mat', 'r')
        elif data_denoise == False and deconv == 'canon':
            f = h5py.File(basedir + 'data/ROI_' + str(r) + str(kfold_data) + 'group_hyp_step2_radius_3_step_2_HOA_cortex_group' + str(curkfold + 1) + '_sub.mat', 'r')
        else:
            f = h5py.File(basedir + 'data/ROI_' + str(r) + str(kfold_data) + 'group_hyp_step2_radius_3_step_2_HOA_cortex_group' + str(curkfold + 1) + '_sub_nr.mat', 'r')

        data = {}
        for k, v in f.items():
            data[k] = np.array(v)
        f.close()
        data = data['x']

        if deconv == 'canon':
            #deconvolve the data
            data = data_deconvolution_canon(data, TR)

        if deconv == 'none' and data_denoise == True:
            data = np.mean(data, axis=2)

    invox = np.nonzero(data[0,:])[0]
    if invox.shape[0] > 0:
        data = data[:, invox]
    else:
        data = np.zeros((np.shape(data)[0], 1))

    return data, hrf_shape



