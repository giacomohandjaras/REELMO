function [results,pvals,beta,null_distro] = maba_measure_Encoding_singlesub(data,model,verbosity,censored_timepoints,statistics,permutation_schema)
%
% [results,pvals,beta,null_distro] = maba_measure_Encoding_singlesub(data,model,verbosity,censored_timepoints,statistics,permutation_schema)
%
% The function uses a 'model' 2D matrix (tps*dims) to predict each column of
% 'matrix', a 2D array of brain activity (voxels*tps). 
% The prediction is performed using OLS multiple linear regression.
% Additional arguments are:
%   a) verbosity: to write debug msg on screen ('yes' | 'no')
%   b) censored_timepoints: an array of tps to censor (e.g., [], if you
%   want to use all the timepoints)
%   c) a string defining the 'statistics' to be used ('none' | 'parametric'
%   | 'permutation')
%   d) a struct 'permutation_schema' containing informations for the
%   permutation test (see create_null_distro_*)
% 
% The function returns:
%   a) 'results' vector (voxels*1) of R^2 coefficients
%   b) 'pvals' vector (voxels*1) of p-values (evaluated with parametric
%   statistics or non-parametric permutation tests)
%   c) 'beta' matrix (voxels*dims) of beta coefficients
%   d) 'null_distro' matrix (voxels*perms) of R^2 coefficients in case of
%   statistics done using the permutation test
%
% authored by giacomo.handjaras@imtlucca.it
%

%% prepare the environment
%%%%% Let's do a backup of the data that will be used in case of a permutation
if strcmp(statistics,'permutation')
    data_ori=data';
end

%%%%% remove censored timepoints and get dimensionality
data=data';
data(censored_timepoints,:)=[];
[tps,voxels]=size(data);

model(censored_timepoints,:)=[];
predictors=size(model,2);

%%%%% Initialize empty arrays
results=[];
pvals=[];
beta=[];
null_distro=[];

%%%%% let's extract MAD. If it is zero means no signal (mask outside fmri data)!
data_mad=mad(data,0,1); 
data_mad_mask=find(data_mad>0);

%%%%% to suppress warnings in case of Nans 
%parfevalOnAll(@warning,0,'off','all');

%% evaluate the effect
%%%%% Initialize the array with the proper dimensions
results=nan(voxels,1);
pvals=nan(voxels,1);
beta=nan(voxels,predictors);

%%%%% Perform multiple regression
selected_data=data(:,data_mad_mask);
beta_data = model\selected_data;

%%%%% calculate R^2
predicted_data = model*beta_data;
ssres = sum((selected_data-predicted_data).^2);
sstot = sum((selected_data-mean(selected_data)).^2);
results(data_mad_mask) = 1-(ssres./sstot);

%%%%% save betas
beta(data_mad_mask,:) = beta_data';

%% 'parametric' statistics
if strcmp(statistics,'parametric')
    %%%%% calculate pvals using betacdf
    pvals(data_mad_mask) = betacdf(results(data_mad_mask), predictors/2, (tps-predictors-1)/2,'upper');
end


%% 'permutation' statistics
if strcmp(statistics,'permutation')
    permutations=size(permutation_schema.permutation_tps,1);
    %%%%% Pick up the first subject of the perm schema
    permutation_schema.permutation_tps=squeeze(permutation_schema.permutation_tps(:,1,:));

    null_distro=nan(voxels,permutations);
    selected_data=data_ori(:,data_mad_mask);
    selected_data_indx=[0:tps:((size(selected_data,2)*tps)-tps)];
    selected_data_indx=repmat(selected_data_indx,tps,1);

    for p=1:permutations
        %if mod(p,100)==0
        %    disp(sprintf('permutation %d out of %d....',p,permutations));
        %end
        current_results=nan(voxels,1);
        current_permutation_schema=repmat(permutation_schema.permutation_tps(p,:)',1,size(selected_data,2));
        current_permutation_schema=current_permutation_schema+selected_data_indx;
        current_selected_data=selected_data(current_permutation_schema);
        current_selected_data(censored_timepoints,:)=[];

        current_beta = model\current_selected_data;
        current_predicted_data = model*current_beta;
        ssres = sum((current_selected_data-current_predicted_data).^2);
        sstot = sum((current_selected_data-mean(current_selected_data)).^2);
        current_results(data_mad_mask) = 1-(ssres./sstot);
        null_distro(:,p)=current_results;
    end

    null_distro_cat=cat(2,results,null_distro);
    null_distro_cat_tie=(tiedrank(null_distro_cat'))';
    pvals(data_mad_mask)=1-((null_distro_cat_tie(data_mad_mask,1)-1)'./max(null_distro_cat_tie(data_mad_mask,:)'));

end

end
