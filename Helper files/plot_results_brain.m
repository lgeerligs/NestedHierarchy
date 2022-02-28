function [vox_vals] = plot_results_brain(fname, vals, inregions, obj, func, basefile,maskvoxin, pvals, pthreshold)
%vals = the values to plot (only for the regions in inregions)

Nregs = size(obj.centers,1);
obj.value=ones(Nregs,1).*NaN;
obj.value(inregions)=vals;

if strcmp(func, 'mean')
    [vox_vals] = map_to_orig_space(obj);
elseif strcmp(func, 'mode')
    [vox_vals] = map_to_orig_space_discrete(obj);
elseif strcmp(func, 'labels')
    [vox_vals, ~] = map_to_orig_space_allvals(obj);
    vox_vals(vox_vals==0)=NaN;
    vox_vals=squeeze(mode(vox_vals,1));
elseif strcmp(func, 'bin')
    [vox_vals] = map_to_orig_space(obj)>0;
end

if exist('pvals', 'var') && exist('pthreshold', 'var')
    if ~isempty(pvals) &&  ~isempty(pthreshold)
        obj.value=ones(Nregs,1).*NaN;
        obj.value(inregions)=pvals;
        [vox_pvals] = map_to_orig_space(obj);
        vox_vals(vox_pvals>pthreshold)=0;
    end
end

vox_vals=vox_vals.*maskvoxin;
hdr=spm_vol(basefile);
hdr.fname=fname;
hdr.dt=[64 0];
spm_write_vol(hdr,vox_vals);

