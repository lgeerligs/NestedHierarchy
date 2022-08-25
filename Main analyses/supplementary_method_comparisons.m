%This script runs the analyses that are shown in the supplementary materials of the manuscript. 
%It is used to compare the reliability of the different methods

clear all

% rootdir='\\cnas.ru.nl\wrkgrp\STD-Donders-DCC-Geerligs\Cambridge_data\';
rootdir = '/home/lingee/wrkgrp/Cambridge_data/';
basedir=[rootdir 'Movie_HMM/'];
toolboxdir=[rootdir 'Toolboxes/'];
resdir=[basedir 'results_fast_event_complete/'];
cd(resdir)

addpath('/home/lingee/wrkgrp/Cambridge_data/Scripts/Movie_HMM/Final code revision/Helper_files')
addpath([toolboxdir '2015_01_25 BCT/'])
addpath(genpath([toolboxdir 'SPM12/']))
addpath(genpath([toolboxdir 'DMLT-master']));
addpath([toolboxdir 'MatlabCentral']);
rmpath(genpath([toolboxdir 'DMLT-master/external/murphy/']))

[ICA_dat, age,gender,mot,subin,rem, CCID, CBUID]=load_subject_info_func(rootdir);

rfiles={basedir 'HAO_mask_cortex.nii'};

%region name
regname='hyp_step2_radius_3_step_2_HOA_cortex';

load([basedir 'voxel_indices' regname '.mat'] );

Ntime=192;
kfold=15;
TR=2.47;

GSdata2states_decon_hrfest=load(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/GSevents_sl_kfold_data2_maxk100_deconv hrf_estsephyp_young_GSBSstates.mat']);
GSdata2states_decon_hrfest.nstates=double(GSdata2states_decon_hrfest.nstates);
GSdata2states_decon_hrfest.binbounds=GSdata2states_decon_hrfest.bounds>0;
GSdata2states_decon_hrfest.binbounds19=GSdata2states_decon_hrfest.bounds<=21&GSdata2states_decon_hrfest.bounds>0;

GSdata2states_decon=load(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/GSevents_sl_kfold_data2_maxk100_deconv canonsephyp_young_GSBSstates.mat']);
GSdata2states_decon.nstates=double(GSdata2states_decon.nstates);
GSdata2states_decon.binbounds=GSdata2states_decon.bounds>0;
GSdata2states_decon.binbounds19=GSdata2states_decon.bounds19>0;

GSdata2states=load(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/GSevents_sl_kfold_data2_maxk100_sephyp_young_GSBSstates.mat']);
GSdata2states.nstates=double(GSdata2states.nstates);
GSdata2states.binbounds=GSdata2states.bounds>0;
GSdata2states.binbounds19=GSdata2states.bounds19>0;

GSdata2bounds=load(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/GSevents_sl_kfold_data2_maxk100_sephyp_young_GSBSbounds.mat']);
GSdata2bounds.nstates=double(GSdata2bounds.nstates);
GSdata2bounds.binbounds=GSdata2bounds.bounds>0;
GSdata2bounds.binbounds19=GSdata2bounds.bounds19>0;

GSdata2bounds_legacy=load(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/GSevents_sl_kfold_data2_maxk100_legacy_GSBS_sephyp_young_GSBSbounds.mat']);
GSdata2bounds_legacy.nstates=double(GSdata2bounds_legacy.nstates);
GSdata2bounds_legacy.binbounds=GSdata2bounds_legacy.bounds>0;
GSdata2bounds_legacy.binbounds19=GSdata2bounds_legacy.bounds<=20&GSdata2bounds_legacy.bounds>0;

load('/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/inregions')
load([resdir 'maskvoxin_searchlight.mat'])

labels = {'bounds FTo', 'bounds FTs', 'states', 'states DC', 'states DE'};

%% show an example of late peak detection and how this is fixed by GSBS states

%supplementary figure 1

[x,y]=find(diff(squeeze(GSdata2bounds.tdists(inregions,1,:))',1,2)>40);
figure; plot(squeeze(GSdata2bounds.tdists(inregions(1237),1,:))); axis([0 100 0 100]); xlabel('Number of states'); ylabel('T-distance')
saveas(gcf,[resdir 'Tdist_GSBS_bounds_k2_f1_r1237.pdf'])

figure; plot(squeeze(GSdata2states.tdists(inregions(1237),1,:))); axis([0 100 0 100]); xlabel('Number of states'); ylabel('T-distance')
saveas(gcf,[resdir 'Tdist_GSBS_states_k2_f1_r1237.pdf'])
%indeed much higher similarity of t-distance curves, higher reliability of
%state boundaries, more reproducable maps. 

%% compare reliability 
kfold=2;
for r=fliplr(inregions')
    for k=kfold:-1:1
        [GSdata2bounds_legacy.median_SL(r,k),~]= median_state_duration(GSdata2bounds_legacy.binbounds(r,k,:));
        [GSdata2bounds.median_SL(r,k),~]= median_state_duration(GSdata2bounds.binbounds(r,k,:));
        [GSdata2states.median_SL(r,k),~]= median_state_duration(GSdata2states.binbounds(r,k,:));
        [GSdata2states_decon.median_SL(r,k),~]= median_state_duration(GSdata2states_decon.binbounds(r,k,:));
        [GSdata2states_decon_hrfest.median_SL(r,k),~]= median_state_duration(GSdata2states_decon_hrfest.binbounds(r,k,:));
    end
end 

median_vox={};
for k=1:2
    [median_vox{k,1}] = plot_results_brain([resdir 'temp.nii'], GSdata2bounds_legacy.median_SL(inregions,k), inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin);
    [median_vox{k,2}] = plot_results_brain([resdir 'temp.nii'], GSdata2bounds.median_SL(inregions,k), inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin);
    [median_vox{k,3}] = plot_results_brain([resdir 'temp.nii'], GSdata2states.median_SL(inregions,k), inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin);
    [median_vox{k,4}] = plot_results_brain([resdir 'temp.nii'], GSdata2states_decon.median_SL(inregions,k), inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin);
    [median_vox{k,5}] = plot_results_brain([resdir 'temp.nii'], GSdata2states_decon_hrfest.median_SL(inregions,k), inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin);
end


rel_legacy = compute_overlap(4, permute(GSdata2bounds_legacy.binbounds(inregions,:,:),[1,3,2]));
rel_bounds = compute_overlap(4, permute(GSdata2bounds.binbounds(inregions,:,:),[1,3,2]));
rel_states = compute_overlap(4, permute(GSdata2states.binbounds(inregions,:,:),[1,3,2]));
rel_statesdecon = compute_overlap(4, permute(GSdata2states_decon.binbounds(inregions,:,:),[1,3,2]));
rel_statesdecon_hrfest = compute_overlap(4, permute(GSdata2states_decon_hrfest.binbounds(inregions,:,:),[1,3,2]));

signrank(squeeze(rel_bounds(:,1,2)), squeeze(rel_states(:,1,2)))

rel19_legacy = compute_overlap(4, permute(GSdata2bounds_legacy.binbounds19(inregions,:,:),[1,3,2]));
rel19_bounds = compute_overlap(4, permute(GSdata2bounds.binbounds19(inregions,:,:),[1,3,2]));
rel19_states = compute_overlap(4, permute(GSdata2states.binbounds19(inregions,:,:),[1,3,2]));
rel19_statesdecon = compute_overlap(4, permute(GSdata2states_decon.binbounds19(inregions,:,:),[1,3,2]));
rel19_statesdecon_hrfest = compute_overlap(4, permute(GSdata2states_decon_hrfest.binbounds19(inregions,:,:),[1,3,2]));

rel = [squeeze(rel_legacy(:,1,2)) squeeze(rel_bounds(:,1,2)) squeeze(rel_states(:,1,2)) squeeze(rel_statesdecon(:,1,2)) squeeze(rel_statesdecon_hrfest(:,1,2))];
mean(rel)
rel19 = [squeeze(rel19_legacy(:,1,2)) squeeze(rel19_bounds(:,1,2)) squeeze(rel19_states(:,1,2)) squeeze(rel19_statesdecon(:,1,2)) squeeze(rel19_statesdecon_hrfest(:,1,2))];
mean(rel19)

barcolors = linspecer(5); 
figure; b=bar(1:5, mean(rel19), 'FaceColor' ,'flat');
for ii=1:size(b.CData,1)
    b.CData(ii,:)=barcolors(ii,:);
end
set(gca, 'Xticklabels',labels);xtickangle(90) 
axis([0 6 0.5 0.6])
hold on; er=errorbar(1:5, mean(rel19), (std(rel19)/sqrt(length(inregions))));
er.Color=[0 0 0];
er.LineStyle='none';
saveas(gcf,[resdir 'compare_methods_overlap_folds_rel19.pdf'])


figure; b=bar(1:5, mean(rel), 'FaceColor' ,'flat');
for ii=1:size(b.CData,1)
    b.CData(ii,:)=barcolors(ii,:);
end
set(gca, 'Xticklabels',labels);xtickangle(90) 
axis([0 6 0.55 0.65])
hold on; er=errorbar(1:5, mean(rel), (std(rel)/sqrt(length(inregions))));
er.Color=[0 0 0];
er.LineStyle='none';
saveas(gcf,[resdir 'compare_methods_overlap_folds_rel.pdf'])


rel_nstates(1)=corr(GSdata2bounds_legacy.median_SL(inregions,1), GSdata2bounds_legacy.median_SL(inregions,2));
rel_nstates(2)=corr(GSdata2bounds.median_SL(inregions,1), GSdata2bounds.median_SL(inregions,2));
rel_nstates(3)=corr(GSdata2states.median_SL(inregions,1), GSdata2states.median_SL(inregions,2));
rel_nstates(4)=corr(GSdata2states_decon.median_SL(inregions,1), GSdata2states_decon.median_SL(inregions,2));
rel_nstates(5)=corr(GSdata2states_decon_hrfest.median_SL(inregions,1), GSdata2states_decon_hrfest.median_SL(inregions,2));

rel_nstates_vox(1)=corr(median_vox{1,1}(median_vox{1,1}>0),median_vox{2,1}(median_vox{2,1}>0), 'rows', 'complete');
rel_nstates_vox(2)=corr(median_vox{1,2}(median_vox{1,2}>0),median_vox{2,2}(median_vox{2,2}>0), 'rows', 'complete');
rel_nstates_vox(3)=corr(median_vox{1,3}(median_vox{1,3}>0),median_vox{2,3}(median_vox{2,3}>0), 'rows', 'complete');
rel_nstates_vox(4)=corr(median_vox{1,4}(median_vox{1,4}>0),median_vox{2,4}(median_vox{2,4}>0), 'rows', 'complete');
rel_nstates_vox(5)=corr(median_vox{1,5}(median_vox{1,5}>0),median_vox{2,5}(median_vox{2,5}>0), 'rows', 'complete');

figure; b=bar(1:5, rel_nstates, 'FaceColor' ,'flat');
for ii=1:size(b.CData,1)
    b.CData(ii,:)=barcolors(ii,:);
end
set(gca, 'Xticklabels',labels);xtickangle(90) 
axis([0 6 0.4 0.65])
saveas(gcf,[resdir 'compare_methods_state_length_sim_folds.pdf'])


%% show that estimating the delays does not make a difference

plot_results_brain([resdir 'estimated_HRF_delays_f1.nii'], squeeze(GSdata2states_decon_hrfest.hrf_shape(inregions,1,2)), inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin);

dat=squeeze(GSdata2states_decon_hrfest.hrf_shape(inregions,1,2));
figure; hist(dat(:),100); xlabel('Estimated HRF peak'); 
ylabel('Frequency'); title('HRF estimation - group 1')
saveas(gcf,[resdir 'HRF_delays1.pdf'])


dat=squeeze(GSdata2states_decon_hrfest.hrf_shape(inregions,1,2));
figure; hist(dat(:),100); xlabel('Estimated HRF peak'); 
ylabel('Frequency'); title('HRF estimation - group 2')
saveas(gcf,[resdir 'HRF_delays2.pdf'])


corr(squeeze(GSdata2states_decon_hrfest.hrf_shape(inregions,1,2)), squeeze(GSdata2states_decon_hrfest.hrf_shape(inregions,2,2)))
%% similarity
all_methods=[GSdata2bounds_legacy.median_SL(inregions,1) GSdata2bounds.median_SL(inregions,1) GSdata2states.median_SL(inregions,1) GSdata2states_decon.median_SL(inregions,1) GSdata2states_decon_hrfest.median_SL(inregions,1)];

method_sim(1)=corr(GSdata2bounds_legacy.median_SL(inregions,1), GSdata2bounds.median_SL(inregions,1));
method_sim(2)=corr(GSdata2bounds.median_SL(inregions,1), GSdata2states.median_SL(inregions,1));
method_sim(3)=corr(GSdata2states.median_SL(inregions,1), GSdata2states_decon.median_SL(inregions,1));
method_sim(4)=corr(GSdata2states_decon.median_SL(inregions,1), GSdata2states_decon_hrfest.median_SL(inregions,1));

corr(GSdata2states_decon.median_SL(inregions,1), GSdata2states_decon_hrfest.median_SL(inregions,1))
corr(GSdata2states_decon.median_SL(inregions,2), GSdata2states_decon_hrfest.median_SL(inregions,2))

corr(median_vox{1,4}(median_vox{1,4}>0), median_vox{1,5}(median_vox{1,5}>0))
corr(median_vox{2,4}(median_vox{2,4}>0), median_vox{2,5}(median_vox{2,5}>0))