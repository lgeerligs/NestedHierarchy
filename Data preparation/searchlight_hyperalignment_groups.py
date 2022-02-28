# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 14:15:55 2018

@author: lg03
"""

import os
import numpy as np
from mvpa2.suite import *
from nibabel import load, save
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt
from scipy import stats
import time

#load subject information
exec(open('load_subject_info_for_statedetection.py').read())
datadir = '/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/Data_for_Donders/'
kfold = 5 #15
subin = np.nonzero(age <= 50)[0]
IDs = CBUID[subin]

start = time.time()
kf = KFold(n_splits=kfold, shuffle=True, random_state=1)

maskname = datadir + 'data_plus_GM_mask.nii'

def compute_ISS(run_datasets):
    subsim = np.empty((len(run_datasets), len(run_datasets)))
    for i in range(len(run_datasets)):
        for j in range(i, len(run_datasets)):
            subsim[i, j] = np.mean(np.mean(np.multiply(run_datasets[i][:,1], run_datasets[j][:,1]), 0))
            subsim[j, i] = subsim[i, j]
    refsub = np.argmax(np.mean(subsim, 0))
    return refsub

count = -1
for train_index, test_index in kf.split(np.arange(0, np.shape(subin)[0])):
    count = count + 1
    datasets = []
    for idx in test_index:
        name = IDs[idx][0] + '_s0w_ME_denoised.nii'
        datasets.append(name)

    run_datasets = []
    for i in datasets:
        print(i)
        fds = fmri_dataset(samples=datadir + i, mask=maskname)
        fds = fds[0:192, :]
        zscore(fds, chunks_attr=None, param_est=None)
        run_datasets.append(fds)

    refsub = compute_ISS(run_datasets)

    hyper = Hyperalignment(level1_equal_weight=True)
    slhyper = SearchlightHyperalignment(radius=3, sparse_radius=2, ref_ds=refsub, compute_recon=False, hyperalignment=hyper)

    slhypmaps = slhyper(run_datasets)
    ds_hyper = [h.forward(sd) for h, sd in zip(slhypmaps, run_datasets)]

    for i in range(len(run_datasets)):
        img = map2nifti(ds_hyper[i])
        save(img, datadir + str(kfold) + 'group_hyperalignment/' + datasets[i][0:9] + '_hyperaligned_g' + str(count) + '.nii')

    print('It took', time.time() - start, 'seconds.')

