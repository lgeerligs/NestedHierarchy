function m = map_to_orig_space_discrete(obj)
% map searchlight results back to original space
% this function is used for data in which overlapping searchlight indices can have different values, report
% the value that occurs most frequently. 

% this code was created by adapting code from the DMLT toolbox https://github.com/distrep/DMLT

% s = obj.spheres;
v = obj.value;
o = obj.original;

nums=max(v);

mt = zeros([nums obj.indims]);

for c=1:length(o)
    
%     if bmask
% %         sidx = msk(s{c});
%     else
%         sidx = s{c};
%     endsize
    sidx = o{c};
    
    if ~isnan(v(c))&&v(c)~=0
        val=v(c);
        mt(val, sidx) = mt(val, sidx) + 1;

    end
    
end

[h,m]=max(mt,[],1);
m(h==0)=0;
m=squeeze(m);




