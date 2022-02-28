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

kfold=1;
Ntime=192;
TR=2.47;

GSdata=load(['/home/lingee/wrkgrp/Cambridge_data/Movie_HMM/results_fast_event_complete/GSevents_sl_kfold_data' num2str(kfold) '_maxk100_deconv canonsephyp_young_GSBSstates.mat']);
GSdata.nstates=double(GSdata.nstates);
GSdata.binbounds=GSdata.bounds>0;

Nregs=length(GSdata.nstates);
load([resdir 'maskvoxin_searchlight.mat'])

%% settings for functional network identification

inname=['GS_boundary_bin_cmat_' num2str(kfold) '_folds'];
gamma=[1:0.1:3];
rep=100;
addpath(genpath([toolboxdir '2019_03_03_BCT']))


%% compute functional connectivity estimates

input_all=zeros([kfold, length(inregions),length(inregions)]);
for k=1:kfold
    input_all(k,:,:)=compute_overlap(1,squeeze(GSdata.binbounds(inregions,k,2:Ntime)));
end

%% run consensus partitioning
parpool(30)

%initialize the variables for the intermediate solutions
temp_partition=zeros(length(gamma),length(inregions), kfold, rep);
ran_partition=zeros(length(gamma),length(inregions), kfold, rep);

%first identify the network partition for each value of gamma, each
%participant group and each repetition. 
for k=1:kfold
    input=squeeze(input_all(k,:,:));
    for n=1:length(gamma)
        g=gamma(n);
        disp([g k])
        parfor r=1:rep
            if sum(ran_partition(n,:,k,r))==0
                [M05,Q05]=community_louvain(input,g,[],'negative_asym');
                temp_partition(n,:,k,r)=M05;
                ran_partition(n,:,k,r)=M05(randperm(length(M05)));
            end
        end
    end
    save([basedir 'results_fast_event_complete/temp_modules_louvain_' inname], 'temp_partition', 'ran_partition');
end

%next run consensus partitioning for each value of gamma
fin_partition=zeros(length(gamma),length(inregions));
fin_partition_noth=zeros(length(gamma),length(inregions));
msim=zeros(length(gamma),1);
msimR=zeros(length(gamma),1);
for n=1:length(gamma)
    g=gamma(n);
    disp(n)
    %compute the agreement matrix over all repetitions of the partitioning
    D=agreement(squeeze(temp_partition(n,:,:)))./(rep*kfold);
    msim(n)=entropy(D(:));
    % the threshold that is used for the agreement matrix is based on permuted labels
    DR=agreement(squeeze(ran_partition(n,:,:)))./(rep*kfold);
    msimR(n)=entropy(DR(:));
    %find a consensus partitioning
    mod = consensus_und(D,mean(DR(:)),rep,g);
    fin_partition(n,:)=mod;
end
save([basedir 'results_fast_event_complete/modules_louvain_gamma_' inname], 'msim', 'msimR', 'fin_partition');

%% compute simularity to Power networks and identify the network sizes and number of networks for each value of gamma

[power] = load_network_labels_power(rootdir,basedir, inregions, obj);
load([basedir 'results_fast_event_complete/modules_louvain_gamma_' inname '.mat'])

powersim=zeros(length(gamma), 1); 
powersimBIG=zeros(length(gamma), 1); 

nummods=zeros(length(gamma), 1);
numbigmods=zeros(length(gamma), 1);
numsingle=zeros(length(gamma), 1);
bignetworks=zeros(size(fin_partition));

for n=1:length(gamma)
    disp(n)
    tempnet=zeros(length(fin_partition(n,:)),1);
    %adapt network labels to only include the large networks
    count=0;
    for i=1:max(fin_partition(n,:))
        if length(find(fin_partition(n,:)==i))>100
            count=count+1;
            tempnet(fin_partition(n,:)==i)=count;
        end
        bignetworks(n,:)=tempnet;
    end
    
    inpower=intersect(find(fin_partition(n,:)>0), find(power>1));
    powersim(n)=AMI(fin_partition(n,inpower),power(inpower)');
    inpower=intersect(find(bignetworks(n,:)>0), find(power>1));
    powersimBIG(n)=AMI(bignetworks(n,inpower),power(inpower)');
    
    nummods(n)=length(unique(fin_partition(n,:)));
    numbigmods(n)=length(unique(bignetworks(n,:)));
    modl=zeros(max(fin_partition(n,:)),1);
    for j=1:max(fin_partition(n,:))
        modl(j)=sum(fin_partition(n,:)==j);
    end
    numsingle(n)=length(find(modl==1));
end


%% select optimal partition

figure; plot(powersim)
hold on; plot(powersimBIG)
legend('power', 'powerBIG')
title('AMI with existing parcellations')
[v,g]=max(powersim);

%project optimal the parcellation to a nifti file (figure 3A)
g=9;
gname=num2str(gamma(g).*10);

networks=fin_partition(g,:);
networks_large=bignetworks(g,:);
plot_results_brain([resdir 'networks_large.nii'], networks_large, inregions, obj, 'labels', [basedir 'GM50mask_allsubs.nii'],maskvoxin)

save([basedir 'results_fast_event_complete/final_modules' gname], 'networks', 'networks_large');  
%% compare to each network to the Power networks and the aDMN and pDMN and provide the data for table S1
[pDMN, aDMN] = load_network_labels_apDMN(rootdir,basedir, inregions, obj);

network_powersim=[];
for i=1:max(networks_large)
    for j=2:max(power)
        val=(networks_large'==i)+(power'==j);
        sim=sum(val>1.1)/sum(networks_large'==i);
        network_powersim(i,j)=sim;
    end
    ind=find(networks_large'==i);
    network_sim_aDMN(i)=mean(aDMN(ind)>0);
    network_sim_pDMN(i)=mean(pDMN(ind)>0);
end

top_sim={}; loc_sim={};
for i=1:max(networks_large)
    [val, loc]=sort(network_powersim(i,:), 'descend');
    if val(1)>2*val(2)
        top_sim{i}=val(1);
        loc_sim{i}=loc(1);
    elseif val(1)>2*val(3)
        top_sim{i}=val(1:2);
        loc_sim{i}=loc(1:2);
    else
        top_sim{i}=val(1:3);
        loc_sim{i}=loc(1:3);
    end
end


top_sim=[];loc_sim=[];
for i=1:max(networks_large)
    [val, loc]=sort(network_powersim(i,:), 'descend');
    top_sim(i,:)=val(1:3);
    loc_sim(i,:)=loc(1:3);
end
