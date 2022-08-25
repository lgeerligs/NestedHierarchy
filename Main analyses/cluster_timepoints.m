clear all

% rootdir='\\cnas.ru.nl\wrkgrp\STD-Donders-DCC-Geerligs\Cambridge_data\';
rootdir = '/home/lingee/wrkgrp/Cambridge_data/';
basedir=[rootdir 'Movie_HMM/'];
toolboxdir=[rootdir 'Toolboxes/'];
resdir=[basedir 'results_fast_event_complete/'];
cd(resdir)

addpath(genpath('/home/lingee/wrkgrp/Cambridge_data/Toolboxes/WSBM_v1.2'))
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
TR=2.47;

load([resdir 'maskvoxin_searchlight.mat'])
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

%load event info
load(['subjective_event_onsets.mat'])
event_onsets_subjective=floor(event_onsets./TR)+1;
event_boundaries_subj=zeros(Ntime,1);
event_boundaries_subj(event_onsets_subjective)=1;
event_boundaries_subj=event_boundaries_subj(2:Ntime);

for kfold=1:2
    GSdata=load([ resdir 'GSevents_sl_kfold_data' num2str(kfold) '_maxk100_deconv canonsephyp_young_GSBSstates.mat']);
    GSdata.nstates=double(GSdata.nstates);
    GSdata.binbounds=GSdata.bounds>0;

    Nregs=length(GSdata.nstates);
        
    %% compute distance

    for kf=1:kfold
        input_overlap=squareform(pdist(squeeze(GSdata.binbounds(inregions,kf,2:Ntime))'));
        input_overlap=max(input_overlap(:))-input_overlap;


        %% WSBM
        LE=[];
        for i=1:10
            [lab{i},model{i}] = wsbm(input_overlap,i,'alpha',0,'parallel',1,'numTrials',1000);
            LE(i,:)=model{i}.LLs;
        end
        figure; plot(LE)

    
        save([resdir 'timepoint_community_grouped_kfold' num2str(kfold) '_fold' num2str(kf) '.mat'], 'LE', 'lab', 'model')

        %% plot ordered connectivity matrix

        %get number of clusters
        [val,i]=max(mean(LE'));

        %order the modules by the average number of boundaries per TR
        data=double(squeeze(GSdata.binbounds(inregions(order_ROIs),kf,2:Ntime)));
        fin_mods=lab{i};
        total_bounds=[];
        for j=1:max(fin_mods)
            total_bounds(j)=mean(mean(data(:,fin_mods==j)));
        end
        [total_bounds_order, order_mods]=sort(total_bounds, 'descend');
        relabel_mods=zeros(size(fin_mods));
        for j=1:max(fin_mods)
            relabel_mods(fin_mods==j)=find(order_mods==j);
        end

        %find the location of the boundaries between modules
        [mods, order]=sort(relabel_mods);
        bounds=[1 ;find(diff(mods)>0)+0.5 ;191];

        %plot average results per network and community
        data=double(squeeze(GSdata.binbounds(inregions(order_ROIs),kf,2:Ntime)));
        data=data(:,order);

        %get averaged data per network
        data_net=zeros(max(ordered_relabeled_nets), max(mods));
        for k=1:max(mods)
            for j=1:max(ordered_relabeled_nets)    
                 %mean number of state per timepoint per network
                data_net(j,k)=squeeze(mean(mean(data(ordered_relabeled_nets==j,mods==k))))/squeeze(mean(mean(data(ordered_relabeled_nets==j,:))));
                
            end
        end
        
        
        %plot the neural states per timepoint
        bounds=find(diff(mods))+0.5;
        bounds=[1; bounds; Ntime];
        ticks=[];
        for j=1:i
            ticks(j)=bounds(j)+((bounds(j+1)-bounds(j))/2);
        end

        %get all data (unaveraged)
        data=double(squeeze(GSdata.binbounds(inregions(order_ROIs),kf,2:Ntime)));
        data(:,event_boundaries_subj==1)=data(:,event_boundaries_subj==1).*2;
        data=data(:,order);

        %get the proportion of event for each module
        number_events=[];
        proportion_events=[];
        for k=1:max(mods)
            [x,y]=find(data(:,mods==k)==2);
            number_events(k)=length(unique(y));
            proportion_events(k)=(length(unique(y))/sum(mods==k))/(length(event_onsets)/Ntime);
        end 
        
        %plot the state transitions per module and network and add the
        %proportion of events
        data_net_events=[data_net; proportion_events];
        figure; imagesc(data_net_events); caxis([0 1]); colormap(hot)
        set(gca, 'Xtick', 1:max(mods), 'Ytick', 1:max(ordered_relabeled_nets)+1, 'Yticklabels', [labels(:); {'events'}])
        print(gcf, '-dpdf', [resdir 'timepoint_community_grouped_event_' num2str(max(mods)) 'kfold' num2str(kfold) '_fold' num2str(kf) '.pdf'])

        % bar plot of neural state and events
        figure; bar(data_net_events)
        set(gca, 'Xtick', 1:max(ordered_relabeled_nets)+1, 'Xticklabels', [labels(:); {'events'}])
        xtickangle(45)
        print(gcf, '-dpdf', [resdir 'timepoint_community_grouped_event_barplot_' num2str(max(mods)) 'kfold' num2str(kfold) '_fold' num2str(kf) '.pdf'])

        
        %plot data averaged per network and per module
        figure; imagesc(data_net);caxis([0 1]); colormap(hot)
        set(gca, 'Xtick', 1:max(mods), 'Ytick', 1:max(ordered_relabeled_nets), 'Yticklabels', labels)
        saveas(gcf, [resdir 'timepoint_community_grouped_' num2str(max(mods)) 'kfold' num2str(kfold) '_fold' num2str(kf) '.pdf'])


        %plot all timepoints and ROIs per module
        fig1 =figure();
        fig1.Renderer='opengl';
        h1=imagesc(data);colormap(hsv)
        data(data==0)=NaN;
        set(h1, 'AlphaData', 1-isnan(data))
        nr=length(order_ROIs);
        hold on;
        for j=1:length(bounds)
            plot(ones(nr,1).*bounds(j), 1:nr, 'k', 'LineWidth', 2);
        end
        for j=1:length(boundaries)
            plot(1:Ntime,ones(Ntime,1).*boundaries(j), 'k--', 'LineWidth', 1);
        end
        set(gca,'Xtick',ticks,'Xticklabels', 1:length(ticks),'Ytick',[])
        xlabel('Community')
        print(gcf, '-dpdf', [resdir 'timepoint_community_' num2str(max(mods)) 'kfold' num2str(kfold) '_fold' num2str(kf) '.pdf'])
    end
end

