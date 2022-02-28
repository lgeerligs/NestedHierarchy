function [centers,spheres,original] = estimate_spheres(obj)
% obtain searchlight spheres
 
%  fields in obj:
%  verbose
%  indims                % mask dimensions
%  exclude = true;       % if false then only the center within a sphere is required to be part of the optional mask

%  radius = 3            % radius of the hypersphere in terms of array elements (diameter will be 1 + 2 * radius).
%  step = 1              % stepsize in terms of array elements    
%     
%  mask                  % optional logical mask of size indims (input features are only those in the mask)
%  neighbours            % a sparse adjacency matrix specifying the neighbourhood structure for irregular data (don't use in conjunction with mask)
%   
%  centers               % centers of each sphere
%  spheres               % the features belonging to each sphere
%  original              % in original space

% this code was created by adapting code from the DMLT toolbox https://github.com/distrep/DMLT

assert(~isempty(obj.indims));

% set default radius
if isempty(obj.radius)
    obj.radius =  max(1,floor(min(obj.indims)./4));
end
rad = obj.radius;

% set default step size
if isempty(obj.step)
    obj.step = max(1,floor(min(obj.indims)./4));
end
stepsize = obj.step;

% identify centers
dd = cell(1,numel(obj.indims));
for c=1:numel(obj.indims)
    dd{c} = 1:stepsize:obj.indims(c);
end

centers = cartprod(dd{:});

% identify centers which are inside the mask
if ~isempty(obj.mask) && ~all(obj.mask(:))
    
    cidx = subv2ind(obj.indims,centers);
    centers = centers(logical(obj.mask(cidx)),:);
    
end

% identify subsets which fall in each sphere

if obj.verbose
    fprintf('estimating %d spheres with radius %g with %g steps for a volume of size [%s]\n',...
        size(centers,1),rad,stepsize,num2str(obj.indims));
end

spheres = cell(size(centers,1),1);

if ~isempty(obj.neighbours)
    
    if obj.verbose
        fprintf('building neighbourhood\n');
    end
    
    % centers as variable indices
    centers = subv2ind(obj.indims,centers);
    
    nidx = cell(size(obj.neighbours,1),1);
    for c=1:length(nidx)
        nidx{c} = find(obj.neighbours(c,:));
    end
    
    for c=1:length(spheres)
        
        if obj.verbose
            fprintf('estimating sphere %d of %d\n',c,size(centers,1));
        end
        
        % explicit neighbourhood structure for irregular data
        % radius becomes path length in the adjacency matrix
        
        spheres{c} = centers(c);
        for j=1:obj.radius
            spheres{c} = unique([spheres{c} cell2mat(nidx(spheres{c})')]);
        end
        
    end
    
else
    
    
    for c=1:length(spheres)
        
        if obj.verbose
            fprintf('estimating sphere %d of %d\n',c,size(centers,1));
        end
        
        % generate potential points (the cube)
        v = cell(1,numel(obj.indims));
        for j=1:numel(obj.indims)
            v{j} = floor(centers(c,j)-rad):ceil(centers(c,j)+rad);
        end
        
        cube = cartprod(v{:});
        
        % reduce to elements inside the cube
        cube = cube(all(cube>0,2) & all(bsxfun(@minus,cube,obj.indims)<=0,2),:);
        
        % reduce to index elements inside the sphere
        spheres{c} = subv2ind(obj.indims,cube(sqrt(sum(bsxfun(@minus,cube,centers(c,:)).^2,2)) <= rad,:));
        
    end
    
    % centers as variable indices
    centers = subv2ind(obj.indims,centers);
    
end

% make spheres fully consistent with mask
if obj.exclude && ~all(obj.mask(:))
    
    midx = find(obj.mask(:));
    
    tmp = spheres;
    for c=1:length(spheres)
        
        % reduce to elements inside the mask
        tmp{c} = spheres{c}(logical(obj.mask(spheres{c})));
        
        % transform to mask indices
        [a,spheres{c}] = ismember(tmp{c},midx);
        
    end
    
    % original sphere indices
    original = tmp;
    
else
    
    original = spheres;
    
end

if obj.verbose
    fprintf('average sphere volume: %g\n',sum(cellfun(@(x)(numel(x)),spheres))/length(spheres));
end
end
