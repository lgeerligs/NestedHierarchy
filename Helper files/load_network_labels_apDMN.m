function [pDMN, aDMN] = load_network_labels_apDMN(rootdir,basedir, inregions, obj)
    
%     if ~exist([basedir 'results_fast_event_complete/searchlight_labels_apDMN.mat'], 'file')

        %create thresholded images with the right dimensions
        spm_imcalc({[basedir 'GM50mask_allsubs.nii'], [rootdir '/templates/Campbell_2013/CKO_yngold_vPCCdPCC_T5_BfMRIbsr_lv2.hdr']}, [rootdir '/templates/Campbell_2013/pDMN.nii'], 'i2>5');
        spm_imcalc({[basedir 'GM50mask_allsubs.nii'], [rootdir '/templates/Campbell_2013/CKO_yngold_MTLPFC_T5_BfMRIbsr_lv2.hdr']}, [rootdir '/templates/Campbell_2013/aDMN.nii'], 'i2<-5');
        hdr1=spm_vol([rootdir '/templates/Campbell_2013/pDMN.nii']);
        hdr2=spm_vol([rootdir '/templates/Campbell_2013/aDMN.nii']);
        
        pDMN=[];
        aDMN=[];
        for r=1:6335
            disp(r)
            %get XYZ coordinates of searchlight
            XYZ = ind2subv(size(obj.mask),obj.original{r})';
           
            [Y1] = spm_get_data(hdr1,XYZ);
            [Y2] = spm_get_data(hdr2,XYZ);
            pDMN(r) = mean(Y1);
            aDMN(r) = mean(Y2);

        end
        pDMN=pDMN(inregions);
        aDMN=aDMN(inregions);

        save([basedir 'results_fast_event_complete/searchlight_labels_apDMN'], 'pDMN', 'aDMN');
%     else
%         load([basedir 'results_fast_event_complete/searchlight_labels_apDMN.mat'], 'pDMN', 'aDMN');
%     end
end