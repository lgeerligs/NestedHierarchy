function m = map_to_orig_space(obj)
% map searchlight results back to original space

% this code was created by adapting code from the DMLT toolbox https://github.com/distrep/DMLT

% s = obj.spheres;
v = obj.value;
o = obj.original;

bmask = ~all(obj.mask(:));
msk = find(obj.mask);

mt = zeros(obj.indims);
n = zeros(obj.indims);
for c=1:length(o)
    
%     if bmask
% %         sidx = msk(s{c});
%     else
%         sidx = s{c};
%     end
    sidx = o{c};
    
    if ~isnan(v(c))
        mt(sidx) = mt(sidx) + v(c);
        n(sidx) = n(sidx) + 1;
    end
    
end

nidx = n(:)~=0;
mt(nidx) = mt(nidx) ./ n(nidx);
mt(n(:)==0) = nan;

m = mt;


