%This script generates all the main figures and statistics reported in the
%paper. 

clear all

% rootdir='\\cnas.ru.nl\wrkgrp\STD-Donders-DCC-Geerligs\Cambridge_data\';
rootdir = '/home/lingee/wrkgrp/Cambridge_data/';
basedir=[rootdir 'Movie_HMM/'];
toolboxdir=[rootdir 'Toolboxes/'];
resdir=[basedir 'results_fast_event_complete/'];
cd(resdir)

addpath(genpath('/home/lingee/wrkgrp/Cambridge_data/Scripts/Movie_HMM/Final code revision'))
addpath([toolboxdir '2015_01_25 BCT/'])
addpath(genpath([toolboxdir 'SPM12/']))
addpath(genpath([toolboxdir 'DMLT-master']));
addpath([toolboxdir 'MatlabCentral']);
rmpath(genpath([toolboxdir 'DMLT-master/external/murphy/']))

[ICA_dat, age,gender,mot,subin,rem, CCID, CBUID]=load_subject_info_func(rootdir);

%load SL voxel indices
load([basedir 'voxel_indiceshyp_step2_radius_3_step_2_HOA_cortex.mat'] );

Ntime=192;
TR=2.47;

GSdata1states_decon=load(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/GSevents_sl_kfold_data1_maxk100_deconv canonsephyp_young_GSBSstates.mat']);
GSdata1states_decon.nstates=double(GSdata1states_decon.nstates);
GSdata1states_decon.binbounds=GSdata1states_decon.bounds>0;

GSdata2states_decon=load(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/GSevents_sl_kfold_data2_maxk100_deconv canonsephyp_young_GSBSstates.mat']);
GSdata2states_decon.nstates=double(GSdata2states_decon.nstates);
GSdata2states_decon.binbounds=GSdata2states_decon.bounds>0;

GSdata15states_decon=load(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/GSevents_sl_kfold_data15_maxk100_deconv canonsephyp_young_GSBSstates.mat']);
GSdata15states_decon.nstates=double(GSdata15states_decon.nstates);
GSdata15states_decon.binbounds=GSdata15states_decon.bounds>0;

%load event info
load(['subjective_event_onsets.mat'])
event_onsets_subjective=floor(event_onsets./TR)+1;
event_boundaries_subj=zeros(Ntime,1);
event_boundaries_subj(event_onsets_subjective)=1;
event_boundaries_subj=event_boundaries_subj(2:Ntime);

%% find included regions and regions with reliability > 0
inregions_temp=find(sum(GSdata15states_decon.nstates,2)>0);

rel_SL15=zeros(length(inregions_temp),15);
rel_p=ones(length(inregions_temp),1).*NaN;

%compute reliability for both methods for each searchlight and compare them
for i=1:length(inregions_temp)
    disp(i)
    for k=1:15
        notk=setdiff(1:15, k);
        rel_SL15(i,k)=corr(squeeze(GSdata15states_decon.binbounds(inregions_temp(i),k,:)), squeeze(mean(GSdata15states_decon.binbounds(inregions_temp(i),notk,:),2)));
    end
    rel_p(i)=signrank(rel_SL15(i,:));
end

[reliability_masked,p_fdr_2] = fdr_bh(rel_p, 0.05, 'dep');

inregions=inregions_temp(reliability_masked);
Nregs=length(GSdata15states_decon.nstates);

save('inregions', 'inregions')

%% load all relevant info about the identified functional networks
load([resdir 'final_modules18.mat'])

%labels and the order of presentation for each of the functional networks
labels={'Aud', 'Vis-early', 'aDMN', 'FPCN', 'Vis-late', 'DAN', 'CON', 'Motor',.... 
    'pDMN', 'sDMN', 'SMN-med', 'SMN-lat'};
order=[8, 11,12,1,2,5,6,7,4,9,10,3];
labels = labels(order);
labels_evs=labels;
labels_evs{13}='all';

relabeled_nets=zeros(size(networks_large));
for i=1:length(labels)
    ind=find(networks_large==order(i));
    relabeled_nets(ind)=i;
end

[ordered_relabeled_nets,order_ROIs]=sort(relabeled_nets);
boundaries=[1 find(diff(relabeled_nets(order_ROIs))) length(relabeled_nets)];

netcolors = linspecer(length(labels));    
newnet_fixorder=zeros(size(networks_large));
for i=1:max(networks_large)
    newnet_fixorder(networks_large==order(i))=i;
end

npcor = max(networks_large);
%% create mean time courses per network 
%exclude the first timepoint, because it can never contain a boundary

%get network representative event timecourse
evs1=zeros(npcor+1, 1, Ntime-1);
evs2=zeros(npcor+1, 2, Ntime-1);
evs15=zeros(npcor+1, 15, Ntime-1);
for i=1:npcor+1
    if i<=npcor
        regs = inregions(find(relabeled_nets==i));
    else
        regs=inregions;
    end
    evs1(i,:)=mean(GSdata1states_decon.binbounds(regs,1,2:Ntime),1);
    for k=1:2
        evs2(i,k,:)=mean(GSdata2states_decon.binbounds(regs,k,2:Ntime),1);
    end
    for k=1:15
        evs15(i,k,:)=mean(GSdata15states_decon.binbounds(regs,k,2:Ntime),1);
    end
end

%% mask image

res=zeros(Nregs,1);
res(inregions_temp)=abs(reliability_masked-1);
obj.value=res;
[invrel_mask] = map_to_orig_space_allvals(obj);

res=zeros(Nregs,1);
res(inregions_temp)=reliability_masked;
obj.value=res;
[coverage_mask] = map_to_orig_space_allvals(obj);

maskvoxin=squeeze(nanmean(coverage_mask,1)>=0.5)&squeeze(nanmean(invrel_mask,1)<0.0000001);
maskvoxout=squeeze(nanmean(coverage_mask,1)<0.5)|squeeze(nanmean(invrel_mask,1)>0);
maskvoxoutrel=squeeze(nanmean(coverage_mask,1)>=0.5)&squeeze(nanmean(invrel_mask,1)>0);

savename_map=[resdir 'coverage.nii'];
hdr=spm_vol([basedir 'GM50mask_allsubs.nii']);
hdr.fname=savename_map;
hdr.dt=[64 0];
spm_write_vol(hdr,[maskvoxin.*2+maskvoxout]);

savename_map=[resdir 'maskvoxoutrel.nii'];
hdr=spm_vol([basedir 'GM50mask_allsubs.nii']);
hdr.fname=savename_map;
hdr.dt=[64 0];
spm_write_vol(hdr,[maskvoxoutrel]);

save([resdir 'maskvoxin_searchlight.mat'],'maskvoxin', 'maskvoxout', 'inregions')

%% plot reliability

plot_results_brain([resdir 'reliability_inregions.nii'], mean(rel_SL15(reliability_masked,:),2), inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'], maskvoxin);
plot_results_brain([resdir 'reliability_inregions_k2.nii'], mean(rel_SL2(reliability_masked,:),2), inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'], maskvoxin);


%% convert number of state to state durations
%results for figure 1, figure S1 

median_SL1=zeros([Nregs, 1]);
median_SL2=zeros([Nregs, 2]);
iqr_SL2=zeros([Nregs, 2]);
% tel(1:2)=1;

for r=inregions'
    disp(r)
    [median_SL1(r),SL_all]= median_state_duration(GSdata1states_decon.binbounds(r,1,:));
    for k=1:2
        [median_SL2(r,k),SL_all]= median_state_duration(GSdata2states_decon.binbounds(r,k,:));
        iqr_SL2(r,k)=iqr(SL_all);
    end
end

%for each fold plot median, iqr
for k=1:2
    [median_voxs{k}] = plot_results_brain([resdir 'median_SL_2fold_f' num2str(k) '.nii'], median_SL2(inregions,k), inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin);
    [iqr_med_voxs{k}] = plot_results_brain([resdir 'iqr_vs_med_SL_2fold_f' num2str(k) '.nii'], iqr_SL2(inregions,k)./median_SL2(inregions,k), inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin);
end


%compute the reliability of state duration across voxels
corr(median_voxs{1}(median_voxs{1}>0), median_voxs{2}(median_voxs{2}>0), 'rows', 'complete')
corr(iqr_med_voxs{1}(median_voxs{1}>0), iqr_med_voxs{2}(median_voxs{2}>0), 'rows', 'complete')

%most frequent values
prctile(median_voxs{1}(median_voxs{1}>0),90)*TR
prctile(median_voxs{1}(median_voxs{1}>0),10)*TR
prctile(median_voxs{2}(median_voxs{2}>0),90)*TR
prctile(median_voxs{2}(median_voxs{2}>0),10)*TR
prctile(median_voxs{1}(median_voxs{1}>0),0)*TR
prctile(median_voxs{1}(median_voxs{1}>0),100)*TR
prctile(median_voxs{2}(median_voxs{2}>0),0)*TR
prctile(median_voxs{2}(median_voxs{2}>0),100)*TR


%compute the reliability of state duration across SLs
corr(median_SL2(inregions,1), median_SL2(inregions,2))
corr(iqr_SL2(inregions,1)./median_SL2(inregions,1), iqr_SL2(inregions,2)./median_SL2(inregions,2))

% is reliability correlated with median state duration?
corr(mean(rel_SL15(reliability_masked,:),2), median_SL1(inregions,:))

%how does reliability - boundary overlap vary with sample size?
rel_SL15_all=compute_overlap(4,permute(GSdata15states_decon.binbounds(inregions_temp,:,:),[1,3,2]));
rel_SL2=compute_overlap(4,permute(GSdata2states_decon.binbounds(inregions_temp,:,:),[1,3,2]));

ind=find(triu(ones(15),1)==1);
mean(rel_SL2(reliability_masked,1,2),1)
mean(mean(rel_SL15_all(reliability_masked,ind),2))



%% make network timescale nifti images for figure 4 
for i=1:max(relabeled_nets)
    plot_results_brain([resdir labels{i} 'network.nii'], ones(sum(relabeled_nets==i),1), inregions(find(relabeled_nets==i)), obj, 'bin', [basedir 'GM50mask_allsubs.nii'],maskvoxin);
    plot_results_brain([resdir 'timescale' labels{i} 'network_f1.nii'], median_SL1(inregions(relabeled_nets==i),1), inregions(find(relabeled_nets==i)), obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin);
end

%% make the correlation matrices shown in figure S3

input_overlap_BG=compute_overlap(5,squeeze(GSdata2states_decon.binbounds(inregions,:,2:Ntime)));
input_overlap=compute_overlap(1,squeeze(GSdata1states_decon.binbounds(inregions,1,2:Ntime)));
input_cor=corr(squeeze(GSdata1states_decon.binbounds(inregions,1,2:Ntime))');
input_ms=corr(squeeze(GSdata1states_decon.meansig(inregions,1,:))');
input_statelen=repmat(median_SL1(inregions,1), [1 length(inregions)])-repmat(median_SL1(inregions,1), [1 length(inregions)])';

plot_imagesc_boundaries([resdir 'overlap_states_BG.pdf'], input_overlap_BG(order_ROIs,order_ROIs), boundaries, [-0.6 0.6])
plot_imagesc_boundaries([resdir 'overlap_states.pdf'], input_overlap(order_ROIs,order_ROIs), boundaries, [-0.6 0.6])
plot_imagesc_boundaries([resdir 'cmat_states.pdf'], input_cor(order_ROIs,order_ROIs), boundaries, [-0.5 0.5])
plot_imagesc_boundaries([resdir 'cmat_meansig.pdf'], input_ms(order_ROIs,order_ROIs), boundaries, [-0.8 0.8])
plot_imagesc_boundaries([resdir 'distmat_nstate.pdf'], abs(input_statelen(order_ROIs,order_ROIs)), boundaries, [-8 8])

%quantify similarity 
mat=ones(length(inregions));
indtriu=find(triu(mat,1)==1);
corr(input_overlap(indtriu),input_overlap_BG(indtriu))
corr(input_overlap(indtriu),input_ms(indtriu))
corr(input_overlap(indtriu),abs(input_statelen(indtriu)))
corr(input_overlap(indtriu),input_cor(indtriu))
corr(input_overlap(indtriu),abs(input_ms(indtriu)))

%% test for significance of overlap between regions across the 15 folds - figure 3B

GSdata=GSdata15states_decon;
overlap_regions=zeros([15, length(inregions),length(inregions)]);
for k=1:15
    overlap_regions(k,:,:)=compute_overlap(1,squeeze(GSdata.binbounds(inregions,k,2:Ntime)));
end
overlap_regions(:,eye(length(inregions))==1)=0;

%which regions show significant overlap?
parpool(10)
pval_overlap=zeros(length(inregions));
parfor i=1:length(inregions)
    disp(i)
    p_temp=zeros(length(inregions),1); 
    for j=1:length(inregions)
        if i<j
            p_temp(j)=signrank(overlap_regions(:,i,j),0, 'tail', 'both');
        end     
    end
    pval_overlap(i,:)=p_temp;
end
mat=ones(length(inregions));
ind2=find(triu(mat,1)==1);
[~,p_fdr_overlap]= fdr_bh(pval_overlap(ind2),0.05);

%plot the overlap
plotoverlap = input_overlap;
plotoverlap(pval_overlap>p_fdr_overlap)=NaN;
%remove all connections in which the direction of effects does not match between k=1 and k=15
plotoverlap(sign(input_overlap)~=sign(squeeze(mean(overlap_regions,1))))=NaN;
plotoverlap(eye(size(plotoverlap))==1)=0;

percsig = sum(plotoverlap(:)>0)/((length(inregions)*length(inregions))-length(inregions))*100;

plot_imagesc_boundaries([resdir 'signboundary_overlap.pdf'], plotoverlap(order_ROIs,order_ROIs), boundaries, [-0.6 0.6])




%% show the network timescales (figure 3B)

figure; H=notBoxPlot(median_SL1(inregions(newnet_fixorder>0),1),newnet_fixorder(newnet_fixorder>0), 'markMedian',true);
set([H.data],...
    'MarkerFaceColor',[1,1,1]*0.8,...
    'markerEdgeColor',[1,1,1]*0.4,...
    'MarkerSize',1)

set([H.mu],'color','k')

for ii=1:length(H)
    set(H(ii).sdPtch,'FaceColor',netcolors(ii,:),...
        'EdgeColor','none')
    
    set(H(ii).semPtch,'FaceColor',netcolors(ii,:),...
        'EdgeColor','none') 
end

box on
set(gca,'Xtick',1:length(labels), 'Xticklabel', labels);
xtickangle(45)
saveas(gcf,[resdir 'networks_boxplot.pdf'])


%% make a figure of 2 specific correlation matrices (inlets figure 1)
%This is done in python with deconvolved data. 
%Here I plot the locations of these SLs on the brain

regs=[2610,2885];%,3820, 4615, 4026];
regnames={'visual', 'mPFC', 'test'};%, 'MFG', 'precuneus', 'sup parietal'};

%make a nifti image of the locations of these searchlights
sls=ones(1,Nregs).*NaN;
sls(regs)=1;
plot_results_brain([resdir 'Locations_of_example_correlation_matrices.nii'], sls, 1:Nregs, obj, 'bin', [basedir 'GM50mask_allsubs.nii'],maskvoxin)


%% compute relative and absolute overlap between neural states and events

kfold=15;
GSdata=GSdata15states_decon;

rel_overlap_events=zeros(kfold, length(inregions));
rel_overlap_events_19=zeros(kfold, length(inregions));
abs_overlap_events=zeros(kfold, length(inregions));
corr_events=zeros(kfold, length(inregions));
for k=1:kfold
    rel_overlap_events(k,:) = compute_overlap(2,squeeze(GSdata.binbounds(inregions,k,2:Ntime)),event_boundaries_subj); 
    rel_overlap_events_19(k,:) = compute_overlap(2,squeeze(double(GSdata.bounds19(inregions,k,2:Ntime)>0)),event_boundaries_subj);
    abs_overlap_events(k,:) = compute_overlap(3,squeeze(double(GSdata.binbounds(inregions,k,2:Ntime)>0)),event_boundaries_subj); 
    corr_events(k,:) = corr(squeeze(GSdata.binbounds(inregions,k,2:Ntime))',event_boundaries_subj); 
end
rel_overlap_events1 = compute_overlap(2,squeeze(GSdata1states_decon.binbounds(inregions,1,2:Ntime)),event_boundaries_subj); 
corr_events1 = corr(squeeze(GSdata1states_decon.binbounds(inregions,1,2:Ntime))',event_boundaries_subj); 
rel_overlap_events1_19 = compute_overlap(2,squeeze(double(GSdata1states_decon.bounds19(inregions,1,2:Ntime)>0)),event_boundaries_subj);
abs_overlap_events1 = compute_overlap(3,squeeze(double(GSdata1states_decon.binbounds(inregions,1,2:Ntime)>0)),event_boundaries_subj);


%which regions show significant overlap with events?
pval_rel_events=zeros(length(inregions),1);
pval_rel_events19=zeros(length(inregions),1);
pval_corr_events=zeros(length(inregions),1);
pval_abs_events=zeros(length(inregions),1);
for i=1:length(inregions)
    pval_rel_events(i) = signrank(rel_overlap_events(:,i),0); 
    pval_corr_events(i) = signrank(corr_events(:,i),0); 
    pval_rel_events19(i) = signrank(rel_overlap_events_19(:,i),0); 
    pval_abs_events(i) = signrank(abs_overlap_events(:,i),0); 
end

overlap1={rel_overlap_events1, rel_overlap_events1_19,abs_overlap_events1,corr_events1};
overlap={rel_overlap_events, rel_overlap_events_19,abs_overlap_events,corr_events};
pvals={pval_rel_events,pval_rel_events19,pval_abs_events,pval_corr_events};
name={'overlap_boundaries_events1_rel.nii', 'overlap_boundaries_events1_rel19.nii',.....
    'overlap_boundaries_events1_abs.nii','corr_boundaries_events1.nii'};

for i=1:4
    [~,p_fdr]= fdr_bh(pvals{i},0.05);
    
    %remove Sls in which the effect is not in the same direction in k1 and k15
    plotevents=overlap1{i};
    plotevents(sign(overlap1{i})~=sign(mean(overlap{i}))')=NaN;
    plotpval=pvals{i};
    plotpval(sign(overlap1{i})~=sign(mean(overlap{i}))')=NaN;
    plot_results_brain([resdir name{i}], plotevents, inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin, plotpval, p_fdr)  
end

%% relate events to binary or strength weighted neural states
%compare 2 predictors: 
%     * regular boundaries
%     * boundary strength

bounds15=GSdata15states_decon.binbounds(inregions,:,2:Ntime);
bounds1=GSdata1states_decon.binbounds(inregions,:,2:Ntime);
str15=GSdata15states_decon.strengths(inregions,:,2:Ntime);
str1=GSdata1states_decon.strengths(inregions,:,2:Ntime);

OA_bin15=zeros(15,length(inregions));
OA_str15=zeros(15,length(inregions));
for k=1:15
    OA_bin15(k,:)=compute_overlap(3,squeeze(bounds15(:,k,:)),event_boundaries_subj);
    OA_str15(k,:)=compute_overlap(3,squeeze(str15(:,k,:)),event_boundaries_subj);
end
OA_bin1=compute_overlap(3,squeeze(bounds1(:,1,:)),event_boundaries_subj);
OA_str1=compute_overlap(3,squeeze(str1(:,1,:)),event_boundaries_subj);

%test if OA increases are significantly different from zero
p_val_str_bin=ones(length(inregions),1).*NaN;

parfor r=1:length(inregions)
    disp(i)
    p_val_str_bin(r)=signrank(OA_str15(:,r)-OA_bin15(:,r),0, 'tail', 'both');
end

[p_masked,p_fdr_str_bin]= fdr_bh(p_val_str_bin,0.05);

%remove associations that are not in the same direction
dep1=OA_str1-OA_bin1;
dep15=mean(OA_str15-OA_bin15)';

plotval=dep1;
plotval(sign(dep1)~=sign(dep15))=NaN;
plotval(dep1<0)=NaN;
plotpval=p_val_str_bin;
plotpval(sign(dep1)~=sign(dep15))=NaN;
plotpval(dep1<0)=NaN;

plot_results_brain([resdir 'OA_str-bin.nii'], plotval, inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin, plotpval, p_fdr_str_bin)

%% look at co-occurence and how it affects the absolute overlap

%compare 3 predictors: 
%     * regular boundaries
%     * boundaries multiplied with overlap within the same network
%     * boundaries multiplied with overlap across the whole brain

bounds15=GSdata15states_decon.binbounds(inregions,:,2:Ntime);
bounds1=GSdata1states_decon.binbounds(inregions,:,2:Ntime);
bound_cooc_WB15=zeros(size(bounds15));
bound_cooc_net15=zeros(size(bounds15));
bound_cooc_WB1=zeros(size(bounds1));
bound_cooc_net1=zeros(size(bounds1));
for i=npcor:-1:1
    regs=find(newnet_fixorder==i);
    bound_cooc_net15(regs,:,:)=bounds15(regs,:,:).*repmat(evs15(i,:,:),[length(regs) 1 1]);
    bound_cooc_net1(regs,:,:)=bounds1(regs,:,:).*repmat(evs1(i,:,:),[length(regs) 1 1]);
end
bound_cooc_WB15(:,:,:)=bounds15(:,:,:).*repmat(evs15(13,:,:),[length(inregions) 1 1]);
bound_cooc_WB1(:,:,:)=bounds1(:,:,:).*repmat(evs1(13,:,:),[length(inregions) 1 1]);

OA_bin15=zeros(15,length(inregions));
OA_net15=zeros(15,length(inregions));
OA_WB15=zeros(15,length(inregions));
for k=1:15
    OA_bin15(k,:)=compute_overlap(3,squeeze(bounds15(:,k,:)),event_boundaries_subj);
    OA_net15(k,:)=compute_overlap(3,squeeze(bound_cooc_net15(:,k,:)),event_boundaries_subj);
    OA_WB15(k,:)=compute_overlap(3,squeeze(bound_cooc_WB15(:,k,:)),event_boundaries_subj);
end
OA_bin1=compute_overlap(3,squeeze(bounds1(:,1,:)),event_boundaries_subj);
OA_net1=compute_overlap(3,squeeze(bound_cooc_net1(:,1,:)),event_boundaries_subj);
OA_WB1=compute_overlap(3,squeeze(bound_cooc_WB1(:,1,:)),event_boundaries_subj);


%test if OA increases are significantly different from zero
p_val_net_bin=ones(length(inregions),1).*NaN;
p_val_WB_bin=ones(length(inregions),1).*NaN;
p_val_WB_net=ones(length(inregions),1).*NaN;
parfor r=1:length(inregions)
    disp(i)
    if ~isnan(OA_net15(1,r))
        p_val_net_bin(r)=signrank(OA_net15(:,r)-OA_bin15(:,r),0, 'tail', 'both');
        p_val_WB_net(r)=signrank(OA_WB15(:,r)-OA_net15(:,r),0, 'tail', 'both');
    end
    p_val_WB_bin(r)=signrank(OA_WB15(:,r)-OA_bin15(:,r),0, 'tail', 'both');
end

dep1={OA_net1-OA_bin1, OA_WB1-OA_bin1, OA_WB1-OA_net1};
dep15={mean(OA_net15-OA_bin15)', mean(OA_WB15-OA_bin15)', mean(OA_WB15-OA_net15)'};
pval={p_val_net_bin,p_val_WB_bin,p_val_WB_net};
name={'OA_net-bin.nii','OA_WB-bin.nii','OA_WB-net.nii'};

for i=1:length(name)
   
    [p_masked,p_fdr]= fdr_bh(pval{i},0.05);
    
    %remove associations that are not in the same direction
    plotval=dep1{i};
    plotval(sign(dep1{i})~=sign(dep15{i}))=NaN;
    plotval(dep1{i}<0&dep15{i}<0)=NaN;
    plotpval=pval{i};
    plotpval(sign(dep1{i})~=sign(dep15{i}))=NaN;
    plotpval(dep1{i}<0&dep15{i}<0)=NaN;
    
    plot_results_brain([resdir name{i}], plotval, inregions, obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin, plotpval, p_fdr)
end

%average within and between networks
for n=npcor:-1:1
    netavg(n)=mean(OA_bin1(relabeled_nets==n));
    netavg_net(n)=mean(OA_net1(relabeled_nets==n));
    netavg_WB(n)=mean(OA_WB1(relabeled_nets==n));
    netavg15(:,n)=mean(OA_bin15(:,relabeled_nets==n),2);
    netavg15_net(:,n)=mean(OA_net15(:,relabeled_nets==n),2);
    netavg15_WB(:,n)=mean(OA_WB15(:,relabeled_nets==n),2);
    p_val_avgnet_bin(n)=signrank(netavg15_net(:,n)-netavg15(:,n),0, 'tail', 'both');
    p_val_avgWB_net(n)=signrank(netavg15_WB(:,n)-netavg15_net(:,n),0, 'tail', 'both');
end
%fdr correction and check if sign matches
[p_masked1,p_fdr]= fdr_bh(p_val_avgnet_bin,0.05);
p_masked1(sign(netavg_net-netavg)~=sign(mean(netavg15_net-netavg15)))=0;
[p_masked2,p_fdr]= fdr_bh(p_val_avgWB_net,0.05);
p_masked2(sign(netavg_WB-netavg_net)~=sign(mean(netavg15_WB-netavg15_net)))=0;

figure; bar([netavg; netavg_net; netavg_WB]')
set(gca,'Xtick',1:npcor, 'Xticklabel', labels_evs(1:12));
xtickangle(45)
saveas(gcf,[resdir 'networks_events_cooc_net_WB.pdf'])

%% plot co-occurence in relation to event boundaries

figure; bar(event_boundaries_subj.*0.91, 'FaceColor', 'r')
hold on; 
ylim([0 0.91])
avgb=squeeze(sum(mean(GSdata1states_decon.binbounds(inregions,:,2:end),2),1));
bar(avgb./length(inregions), 'FaceColor', 'b')
saveas(gcf,[resdir 'Boundary_cooc_bar.pdf'])


%% look at interactions on the regional level (searchlight pairs)

parpool(30)
pctRunOnAll warning('off', 'all')

overlap_allregs_both_k1=zeros(length(inregions));
overlap_allregs_k1_i=zeros(length(inregions));
bounds=squeeze(GSdata1states_decon.binbounds(inregions(order_ROIs),:,2:end));
parfor i=1:length(inregions)
    disp(i)
    tic
    toverlap_allregs_both=zeros(length(inregions),1);
    toverlap_allregs_i=zeros(length(inregions),1);
    for j=1:length(inregions)
       toverlap_allregs_both(j)=compute_overlap(3,bounds(i,:)==1&bounds(j,:)==1,event_boundaries_subj);
       toverlap_allregs_i(j)=compute_overlap(3,bounds(i,:)==1&bounds(j,:)==0,event_boundaries_subj);
    end
    overlap_allregs_both_k1(i,:)=toverlap_allregs_both;
    overlap_allregs_k1_i(i,:)=toverlap_allregs_i;
    toc
end

overlap_allregs_both_k15=zeros(15,length(inregions),length(inregions));
overlap_allregs_k15_i=zeros(15,length(inregions),length(inregions));
for k=1:15
    bounds=squeeze(GSdata15states_decon.binbounds(inregions(order_ROIs),k,2:end));
    parfor i=1:length(inregions)
        disp([k i])
        tic
        toverlap_allregs_both=zeros(length(inregions),1);
        toverlap_allregs_i=zeros(length(inregions),1);
        for j=1:length(inregions)
            toverlap_allregs_both(j)=compute_overlap(3,bounds(i,:)==1&bounds(j,:)==1,event_boundaries_subj)
            toverlap_allregs_i(j)=compute_overlap(3,bounds(i,:)==1&bounds(j,:)==0,event_boundaries_subj)
        end
        overlap_allregs_int_k15(k,i,:)=toverlap_allregs_both;
        overlap_allregs_k15_i(k,i,:)=toverlap_allregs_i;
        toc
    end
end

maxoverlap=max(cat(4,overlap_allregs_k15_i, permute(overlap_allregs_k15_i,[1 3 2])),[],4);
overlap_dif_max=overlap_allregs_int_k15-maxoverlap;
clear maxoverlap overlap_allregs_k15_i overlap_allregs_int_k15

%which regions show significant interactions?
pval_both_vs_max=zeros(length(inregions));
parfor i=1:length(inregions)
    disp(i)
    p_temp=zeros(length(inregions),1); 
    for j=1:length(inregions)
        if i<j
            p_temp(j)=signrank(overlap_dif_max(:,i,j),0, 'tail', 'right');
        end     
    end
    pval_both_vs_max(i,:)=p_temp;
end

mat=ones(length(inregions));
ind2=find(triu(mat,1)==1);
[p_masked,p_fdr_overlap_max]= fdr_bh(pval_both_vs_max(ind2),0.05);

maxoverlapk1=max(cat(3,overlap_allregs_k1_i, overlap_allregs_k1_i'),[],3);
plotval_max=overlap_allregs_both_k1-maxoverlapk1;
sig_max=sign(squeeze(mean(overlap_dif_max)))==sign(plotval_max);
sig_max=sig_max.*((pval_both_vs_max+pval_both_vs_max')<p_fdr_overlap_max);
plotval_max=plotval_max.*sig_max;
plotval_max(plotval_max==0)=NaN;

plot_imagesc_boundaries([resdir 'event_overlap_max.pdf'], plotval_max, [], [-0.4 0.4], [],[],[],'k')
plot_results_brain([resdir 'event_overlap_1dim.nii'], mean(plotval_max>0), inregions(order_ROIs), obj, 'mean', [basedir 'GM50mask_allsubs.nii'],maskvoxin)

for i=1:npcor
    for j=1:npcor
        avgval(i,j)=mean(mean(plotval_max(ordered_relabeled_nets==i,ordered_relabeled_nets==j)>0));
    end
end
    
plot_imagesc_boundaries([resdir 'event_overlap_max_network.pdf'], avgval, [], [-0.3 0.3])

perc_sig=sum(plotval_max(:)>0)/(numel(plotval_max)-length(inregions))*100;

indices=find(triu(eye(length(inregions))==0));
mean_effect_size_k1=nanmean(overlap_allregs_both_k1(indices)-maxoverlapk1(indices));
mean_effect_size_k15=nanmean(overlap_dif_max(:,indices),2);

indices=find(plotval_max(:)>0);
sig_effect_size_k1=nanmean(overlap_allregs_both_k1(indices)-maxoverlapk1(indices));
sig_effect_size_k15=nanmean(overlap_dif_max(:,indices),2);

figure; plot(sort(mean_effect_size_k15), 'LineWidth', 2); hold on; yline(mean_effect_size_k1, 'r', 'LineWidth', 2);
saveas(gcf, [resdir 'mean_effect_size_k15_k1.pdf'])
figure; plot(sort(sig_effect_size_k15), 'LineWidth', 2); hold on; yline(sig_effect_size_k1, 'r', 'LineWidth', 2);
saveas(gcf, [resdir 'sig_effect_size_k15_k1.pdf'])

%% check for a significant association between head motion and neural states
kfold=15;
GSdata=GSdata15states_decon; 

folds = load([resdir 'Kfolds' num2str(kfold) '.mat']);
avgmot=zeros(kfold, 192);

%average head motion within folds for each TR
for k=1:kfold
    inds=folds.subin(find(folds.folds==k));
    avgmot(k,:) = mean(mot.rel_rms(inds,1:192),1);
end

%compute similarity of motion timecourses to neural state boundaries
simns = zeros(Nregs,kfold);
for k=1:kfold
    simns(:,k) = corr(squeeze(GSdata.binbounds(:,k,2:end))',avgmot(k,2:end)');
end

%test if correlations are sign different from zero
mot_p=ones(Nregs,1).*NaN;
for r=inregions'
    mot_p(r)=signrank(simns(r,:));
end
[p_masked,p_fdr_2] = fdr_bh(mot_p(inregions), 0.05);

sum(p_masked)

