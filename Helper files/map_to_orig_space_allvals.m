function [m, nums]= map_to_orig_space_allvals(obj)
% map searchlight results back to original space
% this function is used for data in which overlapping searchlight indices can have different values, report
% all values that occur in the overlapping searchlights 

% this code was created by adapting code from the DMLT toolbox https://github.com/distrep/DMLT

% s = obj.spheres;
v = obj.value;
o = obj.original;

bmask = ~all(obj.mask(:));
msk = find(obj.mask);

nums = zeros(obj.indims);

for c=1:length(o)
    sidx = o{c};    
    nums(sidx) = nums(sidx) + 1;
end

mt = ones([max(nums(:)), obj.indims]).*NaN;
n = ones(obj.indims);

for c=1:length(o)
    
%     if bmask
% %         sidx = msk(s{c});
%     else
%         sidx = s{c};
%     end
    sidx = o{c};
    
    if ~isnan(v(c))
        for vox=sidx'
            mt(n(vox),vox) = v(c);
            n(vox) = n(vox) + 1;
        end
    end
end

m = mt;


