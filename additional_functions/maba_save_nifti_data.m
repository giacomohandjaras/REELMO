function [] = maba_save_nifti_data(nifti_file, data, data_labels, metadata, verbosity)
%
% [] = maba_save_nifti_data(nifti_file, data, data_labels, metadata, verbosity)
%
% Given a 'nifti_file' string filename, a 2D 'data' matrix voxels*dims, 
% a string 'data_labels' (a char with dim names separated by spaces), 
% a 'metadata' struct with information related to a nii volume (see
% maba_load_nifti_mask), the function save data in a 3D/4D nifti file.
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


RESULTS_3D_dims=size(data,2);
RESULTS_3D=zeros(metadata.dims(1),metadata.dims(2),metadata.dims(3),RESULTS_3D_dims);
RESULTS_3D=maba_fill_nifti_data(data,RESULTS_3D,metadata.coordinates);

RESULTS_3D_nii=make_nii(RESULTS_3D, metadata.voxel_size, [0 0 0]);
if verbose; fprintf('Saving nifti file: [\b %s]\b\n', nifti_file); end
save_nii(RESULTS_3D_nii, nifti_file);

fprintf('\n[\bIf you are using AFNI/FSL, please fix the nifti header using AFNI]\b\n');
fprintf('3drefit -newid -view tlrc -space MNI -duporigin %s %s\n',metadata.filename, nifti_file);

if iscell(data_labels)
    data_labels=cell2mat(data_labels);
end

if isstring(data_labels) || ischar(data_labels)
    if strlength(data_labels)>0
        fprintf('3drefit -relabel_all_str ''%s'' %s\n\n', data_labels, nifti_file);
    end
end

end