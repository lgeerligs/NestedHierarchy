import numpy as np
from joblib import Parallel, delayed
import scipy.io as io
from run_GSBS import run_state_detection_fixgroup
import pickle
from sklearn.model_selection import KFold

basedir = '/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/'

#load subject information
exec(open('../Helper files/load_subject_info_for_statedetection.py').read())

groupname = 'sephyp_young'
savedir = basedir + 'results_fast_event/'
savedir_final = basedir + 'results_fast_event_complete/'
nregs=6335
maxk=100

subin = np.nonzero(age <= 50)[0]
temp = io.loadmat(basedir + 'results_complete/inregions.mat')
inregions = np.arange(1,nregs)
statewise = False
kfold_data = 2
TR=2.47
deconv = 'none'#'none' #hrf_est; 'none'; canon
legacy_GSBS=True
data_denoise = False

#this is only used to compare the reliability to the original method
if legacy_GSBS:
    methodname = 'legacy_GSBS_'
    statewise = False
    data_denoise = False
    deconv = 'none'
else:
    methodname = ''

if deconv == 'none':
    deconvname = ''
else:
    deconvname = 'deconv ' + deconv

if data_denoise:
    dataname = 'data_denoise'
else:
    dataname = ''

if statewise==False:
    fname = '_GSBSbounds.mat'
elif statewise==True:
    fname = '_GSBSstates.mat'


savename = methodname + dataname + deconvname + groupname + fname

#Parallel(n_jobs=60)(delayed(run_data_deconvolution)(r, kfold_data, basedir, TR) for r in inregions)

Parallel(n_jobs=60)(delayed(run_state_detection_fixgroup)(r, kfold_data=kfold_data, maxk = maxk, savename=savename, savedir=savedir,basedir=basedir, statewise=statewise, overwrite=True, deconv=deconv, TR=TR, data_denoise=data_denoise, legacy_GSBS=legacy_GSBS) for r in inregions)

# collect results
ntime=192
nstates = np.full([nregs,kfold_data], 0).astype(int)
tdists = np.full([nregs,kfold_data, maxk + 2], 0).astype(float)
bounds = np.full([nregs,kfold_data, ntime],0).astype(int)
strengths = np.full([nregs,kfold_data, ntime],0).astype(float)
meansig = np.full([nregs,kfold_data, ntime],0).astype(float)
bounds19 = np.full([nregs,kfold_data, ntime],0).astype(int)
strengths19 = np.full([nregs,kfold_data, ntime],0).astype(float)
hrf_shape = np.full([nregs, kfold_data, 3], 0).astype(float)

# load data
for r in inregions:
    print(r)
    try:
        filename = savedir + 'GSevents_sl' + str(r) + '_kfold_data' + str(kfold_data) + '_maxk'  + str(maxk) + '_' + savename + '.pickle'
        file = open(filename, 'rb')
        res = pickle.load(file)
        nstates[r - 1,:] = res['nstates']
        tdists[r - 1,:,:] = res['tdists']
        bounds[r - 1, :, :] = res['bounds']
        strengths[r - 1, :,:] = res['strengths']
        bounds19[r - 1, :, :] = res['bounds19']
        strengths19[r - 1, :, :] = res['strengths19']
        meansig[r - 1, :,:] = res['meansig']
        hrf_shape[r - 1, :, :] = res['hrf_shape']
    except(IOError, EOFError):
        continue

io.savemat(savedir_final + 'GSevents_sl' + '_kfold_data' + str(kfold_data) + '_maxk' + str(maxk) + '_' + savename,
           {'nstates': nstates, 'tdists': tdists, 'bounds': bounds, 'strengths': strengths, 'meansig':meansig, 'bounds19':bounds19, 'strengths19':strengths19, 'hrf_shape':hrf_shape})

#save .mat file of subjects assigned to folds
kfold=5
folds=np.zeros(np.shape(subin))
kf = KFold(n_splits=kfold, shuffle=True, random_state=1)
count=0
for train_index, test_index in kf.split(np.arange(0, np.shape(subin)[0])):
    count=count+1
    folds[test_index]=count

io.savemat(savedir_final + 'Kfolds' + str(kfold) + '.mat',{'folds': folds, 'subin': subin+1})
