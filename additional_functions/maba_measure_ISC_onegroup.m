function [results,pvals,null_distro] = maba_measure_ISC_onegroup(data,verbosity,tollerance,censored_timepoints,statistics,permutation_schema)
%
% [results,pvals,null_distro] = maba_measure_ISC_onegroup(data,verbosity,tollerance,censored_timepoints,statistics,permutation_schema)
%
% The function perform Inter-Subject Correlation (ISC) in one group of
% participants. 'data' is a 3D array of brain activity (subjects*voxels*tps). 
% ISC is performed using Pearson's correlation coefficient.
% Additional arguments are:
%   a) tollerance: number in % below which a voxel is discarded from analyses 
%   (e.g., not all the voxels are acquired in all participants)
%   b) verbosity: to write debug msg on screen ('yes' | 'no')
%   c) censored_timepoints: an array of tps to censor (e.g., [], if you
%   want to use all the timepoints)
%   d) a string defining the 'statistics' to be used ('none' | 'parametric'
%   | 'non-parametric' | 'permutation')
%   e) a struct 'permutation_schema' containing informations for the
%   permutation test (see create_null_distro_*)
% 
% The function returns:
%   a) 'results' 2D matrix (voxels*pairings_of_subjects) of r coefficients
%   b) 'pvals' vector (voxels*1) of p-values (evaluated with parametric
%   statistics or non-parametric or permutation tests)
%   c) 'null_distro' [MISSING]
%
% authored by giacomo.handjaras@imtlucca.it
%

%% prepare the environment
data(:,:,censored_timepoints)=[];
[subjects,voxels,tps]=size(data);
pairings=subjects*(subjects-1)/2;
tollerance_minimum=floor(subjects*(tollerance/100));

results=[];
pvals=[];
null_distro=[];

%%%%% let's extract MAD. If it is zero means no signal (mask outside fmri data)!
data_mad=mad(data,0,3); 
data_mad_mask=data_mad>0;
data_mad_mask=sum(data_mad_mask,1);

%%%%% to suppress warnings in case of Nans 
parfevalOnAll(@warning,0,'off','all');

%% evaluate the effect
results=nan(voxels,pairings);
results_temp=cell(voxels,1);

parfor v=1:voxels
    if(data_mad_mask(v)>=tollerance_minimum)
        results_temp{v}=1-pdist(squeeze(data(:,v,:)),'correlation');
    end
end

for v=1:voxels
    if(numel(results_temp{v})>0)
        results(v,:)=results_temp{v};
    end
end

pvals=nan(voxels,1);

%% 'parametric' statistics: r-to-z Fisher transform, ttest vs 0, right tail only
if strcmp(statistics,'parametric')
    for v=1:voxels
        if(numel(results_temp{v})>0)
            [~,pvals(v,1)]=ttest(atanh(results(v,:)),0,'tail','right');
        end
    end
end

%% 'non-parametric' statistics: signrank vs 0, right tail only
if strcmp(statistics,'non-parametric')
    for v=1:voxels
        if(numel(results_temp{v})>0)
            [pvals(v,1),~]=signrank(results(v,:),0,'tail','right');
        end
    end
end

%% 'permutation' statistics
if strcmp(statistics,'permutation')
    error('ERROR: PERMUTATION TEST NOT YET SUPPORTED!');
end

end
