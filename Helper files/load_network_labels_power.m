function [power] = load_network_labels_power(rootdir,basedir, inregions, obj)
    
%     if ~exist([basedir 'results_fast_event_complete/searchlight_labels.mat'], 'file')

        %relate each final module to different predefined functional networks
        hdr2=spm_vol([rootdir '/templates/Power - 2011/Power_consensus_dim.nii']);
        addpath(genpath([rootdir '/Toolboxes/DMLT-master']));
        power=[];
        for r=1:6335
            disp(r)
            %get XYZ coordinates of searchlight
            XYZ = ind2subv(size(obj.mask),obj.original{r})';
           
            [Y1] = spm_get_data(hdr2,XYZ);
            if sum(Y1)==0
                power(r)=0;
            else
                Y1=Y1(Y1>1);
                power(r)=mode(round(Y1));
            end

        end
        power=power(inregions);

        save([basedir 'results_fast_event_complete/searchlight_labels'], 'power');
%     else
%         load([basedir 'results_fast_event_complete/searchlight_labels.mat'], 'power');
%     end
end