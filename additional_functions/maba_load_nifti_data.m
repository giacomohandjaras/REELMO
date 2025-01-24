function [volume] = maba_load_nifti_data(nifti_file,verbosity)
%
% [volume] = maba_load_nifti_data(nifti_file,verbosity)
%
% Given a 'nifti_file' string filename referring to a 3D/4D nifti file 
% (e.g., a .nii fMRI run), the function return the 3D 'volume' matrix.
% If you want to open a 3D file of a binary mask (or an atlas), please 
% refer to maba_load_nifti_mask instead.
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


fprintf('Open nifti data [\b %s]\b\n', nifti_file);
data=load_nii(nifti_file);

dimensions= data.hdr.dime.dim(1);
if verbose; fprintf('Dimensions found %d\n', dimensions); end
current_dimensions=data.hdr.dime.dim(2:dimensions+1);

if sum(current_dimensions==1)>0
    fprintf('[\bA singleton dimension was found...this should not happen! ]\b\n');
end

volume=squeeze(data.img);
voxel_size=data.hdr.dime.pixdim(2:4);
if verbose; fprintf('Voxel size %d x %d x %d mm\n', voxel_size(1), voxel_size(2), voxel_size(3)); end

x_size=size(volume,1);
y_size=size(volume,2);
z_size=size(volume,3);
t_size=size(volume,4);

if verbose; fprintf('Matrix size %d x %d x %d x %d \n', x_size, y_size, z_size,t_size); end

end

