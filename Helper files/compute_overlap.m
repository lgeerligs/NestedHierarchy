function overlap = compute_overlap(mode, varargin)
%This function computes the boundary overlap between brain regions or between regions and events
%mode 1 = relative overlap between regions
%mode 2 = relative overlap between regions and events
%mode 3 = absolute overlap between regions and events (scaled wrt the number of boundaries in a region)
%mode 4 = percentage of overlapping boundaries between folds
%mode 5 = relative overlap between regions for two independent folds

%argument 2 is the data per region (region x time) or (region x time x fold if mode ==4)
%argument 3 is the event timecourse (time x 1)

data = varargin{1};

Ntime = size(data,2);
Nregs = size(data,1);
sumbounds = squeeze(sum(data,2));
ns = sumbounds/Ntime;

if mode==1
    expected_overlap = ns*ns'*Ntime;
    maximum_overlap = squeeze(min(cat(3, repmat(sumbounds,[1, Nregs])', repmat(sumbounds,[1, Nregs])),[],3));
    real_overlap = data*data';
    overlap = (real_overlap-expected_overlap)./(maximum_overlap-expected_overlap);
    
elseif mode == 2
    data2 = varargin{2};
    ns2 = ones(Nregs,1).*sum(data2)/(Ntime);
    expected_overlap = ns.*ns2.*Ntime;
    maximum_overlap = squeeze(min(cat(2, sumbounds, repmat(sum(data2),[Nregs,1])),[],2));
    real_overlap = data*data2;
    overlap = (real_overlap-expected_overlap)./(maximum_overlap-expected_overlap);
    
elseif mode == 3
    data2 = varargin{2};
    ns2 = ones(Nregs,1).*sum(data2)/(Ntime);
    expected_overlap = ns.*ns2.*Ntime;
    maximum_overlap = sumbounds;
    real_overlap = data*data2;
    overlap = (real_overlap-expected_overlap)./(maximum_overlap-expected_overlap);
    
elseif mode ==4
    folds = size(data,3);
    overlap = zeros(Nregs,folds, folds);
    for i=1:folds
        for j=1:folds
            if i>j
                real_overlap = sum(squeeze(data(:,:,i)).*squeeze(data(:,:,j)),2);
                overlap(:,i,j) = real_overlap./((sumbounds(:,i)+sumbounds(:,j))/2);
                overlap(:,j,i)=overlap(:,i,j);
            end
        end
    end
    % overlap between two folds (inter-group overlap)
elseif mode == 5
    Ntime=size(data,3);
    data1 = squeeze(data(:,1,:));
    data2 = squeeze(data(:,2,:));
    sumbounds1 = squeeze(sum(data1,2));
    sumbounds2 = squeeze(sum(data2,2));
    ns1 = sumbounds1/Ntime;
    ns2 = sumbounds2/Ntime;
    expected_overlap = ns1*ns2'.*Ntime;
    maximum_overlap = min(repmat(sumbounds1,[1 Nregs]), repmat(sumbounds2,[1 Nregs])');
    real_overlap = data1*data2';
    overlap = (real_overlap-expected_overlap)./(maximum_overlap-expected_overlap);
    overlap = (overlap+overlap')./2;
end
end