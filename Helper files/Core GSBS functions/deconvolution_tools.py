import numpy as np
import os
from scipy.io import loadmat, savemat
from scipy import signal
import h5py

def data_deconvolution_canon(data, TR):
#this function performs deconvolution based on the canonical HRF
    hrf = spm_hrf(TR)
    data_decon_canon = np.zeros(data.shape)

    if data.ndim == 3:
        for i in range(0, data.shape[2]):
            data_decon_canon[:, :, i] = deconv(data[:, :, i], hrf)
        data_decon_canon = np.mean(data_decon_canon, axis=2)
    elif data.ndim == 2:
        data_decon_canon = deconv(data, hrf)
    return data_decon_canon
	
def run_data_deconvolution(r, kfold_data, basedir, TR):
	#this function performs deconvolution after estimating the HRF shape
    #data shape is time*voxel*subject

    for kfold in range(0, kfold_data):
        print(str(r) + str(kfold))
        data_nr = {}
        f = h5py.File( basedir + 'data/ROI_' + str(r) + str(kfold_data) + 'group_hyp_step2_radius_3_step_2_HOA_cortex_group' + str(kfold+1) + '_sub_nr.mat', 'r')
        for k, v in f.items():
            data_nr[k] = np.array(v)
        f.close()
        data_nr = data_nr['x']
        #invox = np.nonzero(np.mean(np.mean(data_nr, axis=0), axis=1))[0]
        invox = np.nonzero(np.mean(data_nr[0,:], axis=2))[0]
        data_nr = data_nr[:,invox,:]

        decon = {}
        if len(invox) > 14:
            data = {}
            f = h5py.File(basedir + 'data/ROI_' + str(r) + str(kfold_data) + 'group_hyp_step2_radius_3_step_2_HOA_cortex_group' + str(kfold + 1) + '_sub.mat', 'r')
            for k, v in f.items():
                data[k] = np.array(v)
            f.close()
            data = data['x']
            data = data[:, invox, :]

            decon['data_denoised'], decon['data_notdenoised'], decon['hrf_shape_all'], decon['hrf_sub'] = data_deconvolution_hrf_estimation(data_nr, data, TR, val=int(str(r) + str(kfold*100)))
            # decon['data'], decon['hrf_shape_all'], decon['hrf_sub'] = data_deconvolution_hrf_estimation(data, data, TR, val=1)

        else:
            decon['data_denoised'] = np.zeros((np.shape(data_nr)[0], np.shape(data_nr)[1]))
            decon['data_notdenoised'] = np.zeros((np.shape(data_nr)[0], np.shape(data_nr)[1]))
            decon['hrf_shape_all'] = np.zeros((3, np.shape(data_nr)[1], np.shape(data_nr)[2]))
            decon['hrf_sub'] = np.zeros((195, np.shape(data_nr)[2]))

        savemat(basedir + 'data/ROI_' + str(r) + str(kfold_data) + 'group_hyp_step2_radius_3_step_2_HOA_cortex_group' + str(kfold+1) + 'hrf_denoised_data.mat', decon)


def data_deconvolution_hrf_estimation(data_denoised, data_original, TR, val=1):
    #compute the hrf parameters with the data in data_est
    #apply the deconvolution to data_decon
    #assume that all voxels within a SL share an HRF

    cdir = os.getcwd() + '/temp/'
    if not os.path.isdir(cdir):
        os.mkdir(cdir)

    #for each voxel
    data_d = np.reshape(data_denoised, (data_denoised.shape[0], data_denoised.shape[1] * data_denoised.shape[2]))
    np.savetxt(cdir + '/data' + str(val) + '.txt', data_d, delimiter=',')
    os.system("rsHRF --ts '" + cdir + "data" + str(val) + ".txt' --estimation canon2dd --output_dir '" + cdir + "' -TR " + str(TR) + " --passband 0.008 0.1 --thr 1 -T 20 --n_jobs 1")
    hrfdata = loadmat(cdir + 'data' + str(val) + '_hrf_deconv.mat')

    #per voxel
    hrf_shape_all = np.reshape(hrfdata['PARA'], (hrfdata['PARA'].shape[0], data_denoised.shape[1], data_denoised.shape[2]))
    hrf_sub = np.reshape(hrfdata['hrfa'], (hrfdata['hrfa'].shape[0], data_denoised.shape[1], data_denoised.shape[2])).mean(axis=1)

    #deconvolve per subject, assuming that each voxel within the searchlight has the same HRF, do this for both the denoised and non-denoised data
    data_decon_denoised = np.zeros(data_denoised.shape)
    data_decon_notdenoised = np.zeros(data_denoised.shape)
    for i in range(0,data_denoised.shape[2]):
        hrf = signal.resample_poly(hrf_sub[:, i], 1, hrfdata['para']['T'][0][0][0][0])
        data_decon_denoised[:,:, i] = deconv(data_denoised[:,:,i], hrf)
        data_decon_notdenoised[:, :, i] = deconv(data_original[:, :, i], hrf)
    data_decon_denoised = np.mean(data_decon_denoised, axis=2)
    data_decon_notdenoised = np.mean(data_decon_notdenoised, axis=2)

    return data_decon_denoised, data_decon_notdenoised, hrf_shape_all, hrf_sub


def deconv(bold_sig_deconv, hrf):
    nobs=bold_sig_deconv.shape[0]
    nvar=bold_sig_deconv.shape[1]
    data_deconv=np.zeros((nobs,nvar))
    for voxel_id in range(nvar):
        H = np.fft.fft(np.append(hrf,np.zeros((nobs - max(hrf.shape), 1))), axis=0)
        M = np.fft.fft(bold_sig_deconv[:, voxel_id])
        data_deconv[:, voxel_id] = np.fft.ifft(H.conj() * M / (H * H.conj() + .1 * np.mean((H * H.conj()))))
    return data_deconv

def get_hrf(xBF):
    dt     = xBF['dt']
    fMRI_T = xBF['T']
    p      = np.array([6, 16, 1, 1, 6, 0, 32], dtype=float)
    p[len(p) - 1] = xBF['len']
    bf = spm_hrf(dt, p, fMRI_T)
    bf = bf[:, np.newaxis]
    return bf

def spm_hrf(RT, P=None, fMRI_T=16):
    from scipy.special import gammaln
    """
    @RT - scan repeat time
    @P  - parameters of the response function (two gamma functions)

    defaults  (seconds)
    %	P[0] - Delay of Response (relative to onset)	    6
    %	P[1] - Delay of Undershoot (relative to onset)     16
    %	P[2] - Dispersion of Response			            1
    %	P[3] - Dispersion of Undershoot			            1
    %	P[4] - Ratio of Response to Undershoot		        6
    %	P[5] - Onset (seconds)				                0
    %	P[6] - Length of Kernel (seconds)		           32

    hrf  - hemodynamic response function
    P    - parameters of the response function
    """
    p = np.array([6, 16, 1, 1, 6, 0, 32], dtype=float)
    if P is not None:
        p[0:len(P)] = P
    _spm_Gpdf = lambda x, h, l: \
        np.exp(h * np.log(l) + (h - 1) * np.log(x) - (l * x) - gammaln(h))
    # modelled hemodynamic response function - {mixture of Gammas}
    dt = RT / float(fMRI_T)
    u = np.arange(0, int(p[6] / dt + 1)) - p[5] / dt
    with np.errstate(divide='ignore'):  # Known division-by-zero
        hrf = _spm_Gpdf(
            u, p[0] / p[2], dt / p[2]
        ) - _spm_Gpdf(
            u, p[1] / p[3], dt / p[3]
        ) / p[4]
    idx = np.arange(0, int((p[6] / RT) + 1)) * fMRI_T
    hrf = hrf[idx]
    hrf = np.nan_to_num(hrf)
    hrf = hrf / np.sum(hrf)
    return hrf