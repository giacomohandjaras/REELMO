function [volume, metadata] = maba_load_nifti_mask(nifti_file,verbosity)
%
% [volume, metadata] = maba_load_nifti_mask(nifti_file,verbosity)
%
% Given a 'nifti_file' string filename referring to a 3D nifti file 
% (e.g., a .nii mask), the function return the 3D 'volume' matrix, and a 
% metadata struct containing dimensions, coordinates of non zero voxels,
% voxel size in mm, etc..
% If you want to open a 4D file, please refer to maba_load_nifti_data
% instead.
% Argument 'verbosity' handles debug msg on screen ('yes' | 'no')
%
% authored by giacomo.handjaras@imtlucca.it
%

verbose=false;

if ~exist('verbosity','var')
	verbosity = '';
end

switch nargin
	case  1
		verbosity = 'yes';    
	case  2	
		mustBeMember(verbosity,{'yes','no'}) 
end
    
if (strcmp(verbosity,'yes')); verbose=true; end

fprintf('Open nifti mask [\b %s]\b\n', nifti_file);
mask=load_nii(nifti_file);

volume=mask.img;
voxel_size=mask.hdr.dime.pixdim(2:4);
if verbose; fprintf('Voxel size %d x %d x %d mm\n', voxel_size(1), voxel_size(2), voxel_size(3)); end

x_size=size(mask.img,1);
y_size=size(mask.img,2);
z_size=size(mask.img,3);
if verbose; fprintf('Matrix size %d x %d x %d\n', x_size, y_size, z_size); end

voxels=0;
coordinates=[];
for x=1:x_size
    for y=1:y_size
        for z=1:z_size
            if mask.img(x,y,z,1)~=0
                voxels=voxels+1;
                coordinates(voxels,:)=[x,y,z,mask.img(x,y,z,:)];
            end
        end
    end
end
if verbose; fprintf('Voxels in the mask %d\n', voxels); end
  
metadata = struct('voxels',voxels,'voxel_size',voxel_size,'dims',[x_size,y_size,z_size],'coordinates',coordinates,'filename',nifti_file);

end

