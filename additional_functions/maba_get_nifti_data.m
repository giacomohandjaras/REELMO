function [data] = maba_get_nifti_data(volume,coordinates)
%
% [data] = maba_get_nifti_data(volume,coordinates)
%
% Given a 3D or 4D fMRI volume and an array of coordinates (2D Matrix:
% voxels*[x,y,z]), the function returns a 2D matrix voxels*tps, 
% where tps is the fourth dimension of the volume.
%
% authored by giacomo.handjaras@imtlucca.it
%

voxels=size(coordinates,1);
tps=size(volume,4);
data=nan(voxels,tps,class(volume));

for v=1:voxels
    data(v,:)=volume(coordinates(v,1),coordinates(v,2),coordinates(v,3),:);
end


end
