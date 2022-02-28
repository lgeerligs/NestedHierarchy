function [obj]=voxel_indices_data_mask(maskimage,dataimage)
%obtains the voxel indices in the data for a given maskimage, even if
%image dimensions are different

%this code is based on the roi_extract function by Rik Henson
%https://github.com/MRC-CBU/riksneurotools/blob/master/Util/roi_extract.m

% read image
VY = spm_vol(dataimage);
[dat,yXYZmm] = spm_read_vols(VY(1)); % Assume all data in same space!
Yinv  = inv(VY(1).mat);

%read mask
VM = spm_vol(maskimage);
Minv = inv(VM.mat);
[YM,mXYZmm] = spm_read_vols(VM);
ROIvals = setdiff(unique(YM(:)),0); 

% Transform ROI XYZ in mm to voxel indices in data:
yXYZind = Yinv(1:3,1:3)*mXYZmm + repmat(Yinv(1:3,4),1,size(mXYZmm,2));
% Transform data XYZ in mm to voxel indices in mask:
mXYZind = Minv(1:3,1:3)*yXYZmm + repmat(Minv(1:3,4),1,size(yXYZmm,2));
% Transform data XYZ in mm to voxel indices in data:
yyXYZind = Yinv(1:3,1:3)*yXYZmm + repmat(Yinv(1:3,4),1,size(yXYZmm,2));

image_data = spm_get_data(VM,mXYZind); 

resized_mask=zeros([size(dat)]);
for r=1:length(ROIvals)
    %get data from the image
    f = find(image_data == ROIvals(r));
    obj.voxels{r}=yyXYZind(:,f);
    obj.original{r}= subv2ind(size(dat),obj.voxels{r}');
    resized_mask(obj.original{r})=r;
end

obj.mask=resized_mask;
