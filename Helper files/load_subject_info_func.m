function [ICA_dat, age,gender,mot,subin,rem, CCID, CBUID]=load_subject_info_func(rootdir)
% load('/home/lingee/wrkgrp/Cambridge_data/FC_analysis/ME_ICA/subs_MEICA.mat')
% load('/home/lingee/wrkgrp/Cambridge_data/FC_analysis/ME_ICA/subs_MEICA_icadata.mat')

    load([rootdir '/FC_analysis/ME_ICA/subs_MEICA.mat'])
    load([rootdir '/FC_analysis/ME_ICA/subs_MEICA_icadata.mat'])


    %select subjects to include
    rem=rejall'./tot';
    subin=find(rem<0.8&rms_tot<(mean(rms_tot)+2*std(rms_tot)));

    ICA_dat.var=var(subin);
    ICA_dat.tot=tot(subin);
    ICA_dat.rejall=rejall(subin);
    ICA_dat.rejnb=rejnb(subin);
    ICA_dat.rem=rem(subin);
    age=age(subin);
    gender=gender(subin);
    mot.rel_rms=rel_rms(subin,:);
    mot.rms_tot=rms_tot(subin);
    mot.rms_max=rms_max(subin);
    mot.rms_skew=rms_skew(subin);
    mot.rms_totlarge=rms_totlarge(subin);
    CCID=CCID(subin);
    CBUID=CBUID(subin);
end