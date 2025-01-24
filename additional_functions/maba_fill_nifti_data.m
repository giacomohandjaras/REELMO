function [volume] = maba_fill_nifti_data(data,volume,coordinates)
%
% [volume] = maba_fill_nifti_data(data,volume,coordinates)
%
% Given a 2D matrix (voxels*tps), a 3D/4D fMRI empty volume, 
% and an array of coordinates (2D Matrix: voxels*[x,y,z]), 
% the function returns the volume filled with data.
%
% authored by giacomo.handjaras@imtlucca.it
%

voxels=size(coordinates,1);
[x,y,z,ndim]=size(data);

for v=1:voxels
    volume(coordinates(v,1),coordinates(v,2),coordinates(v,3),:)=data(v,:);
end


end
