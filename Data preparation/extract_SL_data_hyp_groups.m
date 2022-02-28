% This script extracts data for each searchlight
% The input files needed in this script are the hyperaligned nifti images

%ind2subv, voxel_indices_data_mask, estimate_spheres

clear all

rootdir = '/home/lingee/wrkgrp/Cambridge_data/';
basedir = '/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/';

toolboxdir = '/home/lingee/wrkgrp/Cambridge_data/Toolboxes/';
confounddir = '/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/Data_for_Donders/global_signals/';

addpath(basedir)
addpath(genpath('/home/lingee/wrkgrp/Cambridge_data/Scripts/Movie_HMM/Final code revision/'))
addpath(genpath([toolboxdir 'SPM12/']))
addpath(genpath([toolboxdir 'DMLT-master']));

[ICA_dat, age,gender,mot,subin,rem, CCID, CBUID]=load_subject_info_func(rootdir);
file='_hyperaligned';

TR=2.47;
Ns=192;
kfold=1;

dataname=[num2str(kfold) 'group_hyp_step2'];
datadir = [basedir 'Data_for_Donders/' num2str(kfold) 'group_hyperalignment/'];
%datadir = ['/huge/Linda/' num2str(kfold) 'group_hyperalignment/'];


if kfold>1
    load(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/Kfolds' num2str(kfold) '.mat']) 
else
    load(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/Kfolds2.mat'], 'subin') 
end
CBUID=CBUID(subin);
nsub=length(subin);
%% get list of participant files

clear hdrs data
%list of all subjects files
if kfold<2
    for s=length(CBUID):-1:1
        disp(s)
        if kfold==1
            hdrs{1}{s}=spm_vol([datadir CBUID{s} file '.nii']);
        elseif kfold==0
            hdrs{1}{s}=spm_vol(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/Data_for_Donders/' CBUID{s} '_s0w_ME_denoised.nii']);
        end
        nuisance = load([confounddir CBUID{s} '.mat']);
        m = nuisance.motion(1:Ns,:);
        CSF = nuisance.x.CSFmean(1:Ns);
        WM = nuisance.x.WMmean(1:Ns);
        dm = [zeros(1,6); diff(m,1,1)];
        dc = [zeros(1,2); diff([CSF WM],1,1)];
        regressors{1}(s,:,:) = [CSF WM m dm dc  CSF.^2 WM.^2 m.^2 dm.^2 dc.^2];
    end
elseif kfold>1
    for k=kfold:-1:1
        inds=find(folds==k);
        for s=length(inds):-1:1
            disp([k s])
            hdrs{k}{s}=spm_vol([datadir CBUID{inds(s)} file '_g' num2str(k-1) '.nii']);
            nuisance = load([confounddir CBUID{s} '.mat']);
            m = nuisance.motion(1:Ns,:);
            CSF = nuisance.x.CSFmean(1:Ns);
            WM = nuisance.x.WMmean(1:Ns);
            dm = [zeros(1,6); diff(m,1,1)];
            dc = [zeros(1,2); diff([CSF WM],1,1)];
            regressors{k}(s,:,:) = [CSF WM m dm dc  CSF.^2 WM.^2 m.^2 dm.^2 dc.^2];
        end
    end
end

exdata=hdrs{1}{1}; 

%% make high pass DCT filter
Kall   = spm_dctmtx(Ns,Ns);
nHP = fix(2*(Ns*TR)/(1/0.008) + 1);
filter   = Kall(:,[1:nHP]);

%% define searchlights

%searchlight settings
rfiles=[basedir 'HAO_mask_cortex.nii'];
rnames='HOA_cortex';
radius=3;
step=2;

%region name
regname=[dataname '_radius_' num2str(radius) '_step_' num2str(step) '_' rnames];
if kfold>0
    if exist([basedir 'voxel_indices' regname  '.mat'] )
        load([basedir 'voxel_indices' regname '.mat'] )
    else
        mask=voxel_indices_data_mask(rfiles,exdata);
        mask.mask(mask.mask>0)=1;
        obj.radius=radius;    obj.step=step;
        obj.mask=mask.mask;   obj.verbose=1; obj.indims=size(obj.mask); obj.neighbours=[];
        obj.exclude = true;
        [obj.centers,obj.spheres,obj.original]=estimate_spheres(obj);
        save([basedir 'voxel_indices' regname] ,'obj')
    end
else
    load([basedir 'voxel_indices1group_hyp_step2_radius_' num2str(radius) '_step_' num2str(step) '_' rnames '.mat'] )
end

nregions=length(obj.original);
disp(nregions)


%% loop over the searchlights to extract data
R  = eye(Ns) - filter*pinv(filter);
if kfold==0; kfold=1; end
for k=1:kfold
    disp(k)
    curhdr = hdrs{k};
    clear imgs
    imgs=zeros(length(curhdr), 61,73,61,192);  
    
    for s=1:length(curhdr)
        disp(s)
        imgs(s,:,:,:,:)=spm_read_vols(curhdr{s}(1:192));
    end
    
    for r=1:nregions
        disp([k r])
        
        %get XYZ coordinates of searchlight
        XYZ = ind2subv(size(obj.mask),obj.original{r})';
         
        %load data for this searchlight
        Y1=zeros(length(curhdr), size(XYZ,2), Ns);
        Y_sub=zeros(length(curhdr), size(XYZ,2), Ns);
        
        for i=1:size(XYZ,2)
            Y1(:,i,:)=imgs(:,XYZ(1,i),XYZ(2,i),XYZ(3,i),:);
        end
        for s=1:length(curhdr)
            dat = R*squeeze(Y1(s,:,:))'; 
            Y_sub(s,:,:)=dat';
        end
        Y=squeeze(mean(Y_sub,1));
        parsave([basedir 'data/ROI_' num2str(r) regname '_group' num2str(k)],Y)
        parsave([basedir 'data/ROI_' num2str(r) regname '_group' num2str(k) '_sub'],Y_sub)
        clear Y Y_sub
        Y_sub_nr=zeros(length(curhdr), size(XYZ,2), Ns);
        for s=1:length(curhdr)
            regs = [squeeze(regressors{k}(s,:,:)) filter];
            R1  = eye(Ns) - regs*pinv(regs);
            dat = R1*squeeze(Y1(s,:,:))'; 
            Y_sub_nr(s,:,:)=dat';
        end

        parsave([basedir 'data/ROI_' num2str(r) regname '_group' num2str(k) '_sub_nr'],Y_sub_nr)
    end
end