
function [vals] = load_SL_vals(input, output, inregions, obj)
    
%     if ~exist(output, 'file')

        %relate each final module to different predefined functional networks
        hdr=spm_vol(input);
 
        addpath(genpath(['/home/lingee/wrkgrp/Cambridge_data/Toolboxes/DMLT-master']));
        vals=[];
        for r=1:6335
            disp(r)
            %get XYZ coordinates of searchlight
            XYZ = ind2subv(size(obj.mask),obj.original{r})';
            [Y1] = spm_get_data(hdr,XYZ);
            if sum(Y1)==0
                vals(r)=NaN;
            else
                Y1=Y1(Y1~=0);
                vals(r)=mean(Y1);
            end

        end
        vals=vals(inregions);


        save(output, 'vals')
%     else
%         load(output, 'vals');
%     end
end