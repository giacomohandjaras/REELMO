%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Script by Giacomo Handjaras                                          %%%%%
%%%%% This script performs encoding of a model using multiple regression   %%%%% 
%%%%% in a fMRI dataset pertaining to a single subject.                    %%%%%
%%%%% It requires a mask of grey matter to isolate voxels of interest,     %%%%%
%%%%% a model (the intercept will be automatically added)                  %%%%%
%%%%% and fMRI data of one subjects.                                       %%%%%
%%%%%                                                                      %%%%%
%%%%% Statistical analyses are performed by using:                         %%%%%
%%%%% 'none' => just reporting the R^2 effect size,                        %%%%%
%%%%% 'parametric' => the pvalue associated to the R^2,                    %%%%%
%%%%% 'permutation' => a non-parametric permutation using a permutation    %%%%%
%%%%%                  schema on fMRI signals, defined by the user.        %%%%%
%%%%%                  P-vals are calculated using Pareto right tail       %%%%%
%%%%%                  approximation.                                      %%%%%
%%%%% Correction for multiple comparisons is handled by using FDR and/or   %%%%%
%%%%% non-parametric permutation using FWE correction.                     %%%%%
%%%%% Group-level analysis are handled by using another script             %%%%%
%%%%% (Encoding _multiplesubs.m) by non-parametric combination of single   %%%%%
%%%%% participants results obtained here.                                  %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

addpath('additional_functions/');
addpath('external_functions/');
addpath('external_functions/NIfTI_tools/');

%%%%% voxels to be analyzed (a binary mask)
filename_mask='../REELMO/jojo-rabbit/fmri/derivatives/group/brain_mask.nii.gz'; 

%%%%% the model we are planning to encode: a CSV file (with header) handled through
%%%%% readtable. The rows are timepoints, the columns are features
filename_model='models/visual_model.csv'; 

%%%%% EPI cleaned and smoothed timeseries of the subject of interest %%%%%
subject_name='all'; %%% the name of the subject (use 'all' for sub-all, that is the name of the average subject)
subject_path='../REELMO/jojo-rabbit/fmri/derivatives/'; %%% the absolute path (in BIDS format)
subject_file='_epi_cleaned_reml_2mni.nii.gz'; %%% the actual fMRI preprocessed files to load (in BIDS format) 

%%%%% files to be saved
filename_save=['../REELMO/jojo-rabbit/fmri/derivatives/group/Results_Encoding_visual']; %%% a nifti .nii and a matlab .mat will be created at the end of the script

%%%%% the duration of each run in volumes. Actually the run_durations is summed to obtain the total duration of the timeserie
run_durations=[ 272 407 402 403 410 390 390 413 ] ; %%%run duration in TRs

%%%%% censored time points
run_censored_timepoints=[];

%%%%% perform statistics! It could be 'none' if we aim to obtain only the effect size, 'parametric' if we want to use the p-val coming from R^2 multiple regression, 'permutation' using a non-parametric permutation test 
statistics='parametric';
statistics_correction='fdr95'; %%%%%'fdr95' or 'fdr01'
statistics_permutation_data=''; %%% permutation schema for statistic='permutation'

%%%%% define the precision of the data
data_precision='single';

%%%%% number of voxel bins (to limit memory use and speedup computations)
%%%%% This script works generally better without parallelization, 
%%%%% since it uses "\" to perform multiple regression which already works 
%%%%% in parallel 
voxels_bins=25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's prepare the environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------\n');
fprintf(' PREPARING ENVIRONMENT          \n');
fprintf('--------------------------------\n');

%%%%% open model
fprintf('Opening the encoding model: [\b %s]\b\n', filename_model);
[encoding_model_table]=readtable(filename_model,'FileType', 'text','Encoding','UTF-8','Delimiter',',');
encoding_model=table2array(encoding_model_table);
% encoding_model(:,1)=[];
%%%%% check if the first column is already the intercept
intercept_flag=0;
encoding_model_mad=mad(encoding_model,0,1); %%% let's extract MAD. If it is zero means no variation (intercept?)
if (sum(encoding_model_mad==0)==0)
    fprintf('Adding the intercept...\n');
    encoding_model=cat(2,ones(size(encoding_model,1),1),encoding_model);
    intercept_flag=1;
end

[encoding_model_rows,encoding_model_columns]=size(encoding_model);
fprintf('The model has %d rows (%d expected)\n', encoding_model_rows,sum(run_durations));
fprintf('The model has %d columns/descriptors (including the intercept)\n', encoding_model_columns);


%%%%% create a vector encoding each run
run_tps=nan(sum(run_durations),1);
run_starting=1;
for r=1:numel(run_durations)
    run_ending=run_starting-1+run_durations(r);
    run_tps(run_starting:run_ending)=r;
    run_starting=run_ending+1;
end
run_tps(run_censored_timepoints)=[]; %%%%% censored tps
clear r run_starting run_ending

%%%%% check statistics
statistics_flag = sum(ismember(statistics,{'none','parametric','permutation'}))>0;
if statistics_flag
    fprintf('Statistics defined: [\b %s]\b\n', statistics);
    subject_random_schema=[];
    subject_permutations=[];
    subject_permuted_subjects=[];
else
    error('Statistics "[\b %s]\b" was not supported!\n', statistics);
end

%%%%% open permutation schema
if strcmp(statistics,'permutation')
    fprintf('Opening permutation schema: [\b %s]\b\n', statistics_permutation_data);
    subject_random_schema=load(statistics_permutation_data);
    subject_permutations=size(subject_random_schema.permutation_tps,1);
    fprintf('Number of permutations: %d\n', subject_permutations);
    subject_permuted_subjects=size(subject_random_schema.permutation_tps,2);
    fprintf('Number of permuted subjects supported: %d\n', subject_permuted_subjects);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's open the fMRI mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------\n');
fprintf(' LOADING MASKS                  \n');
fprintf('--------------------------------\n');

[mask_volume,mask_metadata]=maba_load_nifti_mask(filename_mask,'yes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's open the subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------\n');
fprintf(' LOADING fMRI                   \n');
fprintf('--------------------------------\n');

subject_filename=strcat(subject_path,'sub-',subject_name,'/func/','sub-',subject_name,subject_file);
[subject_volume]=maba_load_nifti_data(subject_filename,'yes');
subject_data=maba_get_nifti_data(subject_volume,mask_metadata.coordinates);
clear subject_volume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's measure Goodness of Fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------\n');
fprintf(' EVALUATE THE MODEL             \n');
fprintf('--------------------------------\n');

%%%%% check precision again
if strcmp(data_precision,'single')
    subject_data=single(subject_data);
    encoding_model=single(encoding_model);
else
    subject_data=double(subject_data);
    encoding_model=double(encoding_model);
end

%%%%% To reduce memory usage, we need to split
%%%%% the voxels in a certain number of bins 

mask_voxels_binsize=floor(mask_metadata.voxels/voxels_bins);
mask_voxels_bins=nan(voxels_bins,2);
starting_point=1;
ending_point=mask_voxels_binsize;
for b=1:voxels_bins
    mask_voxels_bins(b,1)=starting_point;
    mask_voxels_bins(b,2)=ending_point;
    starting_point=mask_voxels_bins(b,2)+1;
    ending_point=ending_point+mask_voxels_binsize;
end
mask_voxels_bins(end,2)=mask_metadata.voxels;
clear starting_point ending_point b

%%%%% Now perform a for loop for each bin
subject_results=nan(mask_metadata.voxels,1,data_precision);
subject_results_pvals=nan(mask_metadata.voxels,1,'double');
subject_results_beta=nan(mask_metadata.voxels,size(encoding_model,2),data_precision);
subject_results_null=[];
if strcmp(statistics,'permutation')
    subject_results_null=nan(mask_metadata.voxels,subject_permutations,data_precision);
end

for b=1:voxels_bins
    tic
    fprintf('bin %d out of %d....\n',b,voxels_bins);
    current_subject_data=subject_data(mask_voxels_bins(b,1):mask_voxels_bins(b,2),:);
    [current_subject_results,current_subject_results_pvals,current_subject_results_beta,current_subject_results_null]=maba_measure_Encoding_singlesub(current_subject_data,encoding_model,'yes',run_censored_timepoints,statistics,subject_random_schema);
    subject_results(mask_voxels_bins(b,1):mask_voxels_bins(b,2),:)=current_subject_results;
    subject_results_pvals(mask_voxels_bins(b,1):mask_voxels_bins(b,2),:)=current_subject_results_pvals;
    subject_results_beta(mask_voxels_bins(b,1):mask_voxels_bins(b,2),:)=current_subject_results_beta;
    if strcmp(statistics,'permutation')
        subject_results_null(mask_voxels_bins(b,1):mask_voxels_bins(b,2),:)=current_subject_results_null;
    end
    clear current_*
    toc
end

clear b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------\n');
fprintf(' ASSEMBLING 3D VOLUME           \n');
fprintf('--------------------------------\n');

%%%%% Prepare the data
if strcmp(statistics,'none')
    results=subject_results;
    results_labels=strcat('R^2');
end

if strcmp(statistics,'parametric')
    subject_results_pvals_log=abs(log10(subject_results_pvals));
    subject_results_pvals_log(isnan(subject_results_pvals_log))=0;
    subject_results_pvals_corr=maba_correct_pvals(subject_results_pvals,'yes',statistics_correction);
    subject_results_pvals_corr_log=abs(log10(subject_results_pvals_corr));
    subject_results_pvals_corr_log(isnan(subject_results_pvals_corr_log))=0;
    results=cat(2,subject_results,subject_results_pvals_log,subject_results_pvals_corr_log);
    results_labels=strcat('R^2',{' '},'log_pval-raw_parametric',{' '},'log_pval-',statistics_correction,'_parametric');
    clear subject_results_pvals_log subject_results_pvals_corr_log;
end

if strcmp(statistics,'permutation')
    subject_results_pvals_log=abs(log10(subject_results_pvals));
    subject_results_pvals_log(isnan(subject_results_pvals_log))=0;
    subject_results_pvals_corr=maba_correct_pvals(subject_results_pvals,'yes',statistics_correction);
    subject_results_pvals_corr_log=abs(log10(subject_results_pvals_corr));
    subject_results_pvals_corr_log(isnan(subject_results_pvals_corr_log))=0;
    subject_results_null_perc50=prctile(subject_results_null,50,2);
    subject_results_null_perc50(isnan(subject_results_null_perc50))=0;
    subject_results_null_perc99=prctile(subject_results_null,99,2);
    subject_results_null_perc99(isnan(subject_results_null_perc99))=0;
    results=cat(2,subject_results,subject_results_pvals_log,subject_results_pvals_corr_log,subject_results_null_perc50,subject_results_null_perc99);
    results_labels=strcat('R^2',{' '},'log_pval-raw_nonparametric',{' '},'log_pval-',statistics_correction,'_nonparametric',{' '},'null_effect_50perc',{' '},'null_effect_99perc');
    clear subject_results_pvals_log subject_results_pvals_corr_log subject_results_null_perc50 subject_results_null_perc99;
end

%%%%% Add beta coeffs, and beta coeff labels
subject_results_beta_proc=subject_results_beta;
subject_results_beta_proc(isnan(subject_results_beta_proc))=0;
results=cat(2,results,subject_results_beta_proc);
clear subject_results_beta_proc

if intercept_flag; results_labels=strcat(results_labels,{' '},'intercept'); end
for i=1:numel(encoding_model_table.Properties.VariableNames)
    results_labels=strcat(results_labels,{' '},strrep(encoding_model_table.Properties.VariableNames{i},' ','_'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's save the results in nifti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------\n');
fprintf(' SAVING ENCODING IN 3D          \n');
fprintf('--------------------------------\n');

filename_save_nii=strcat(filename_save,'.nii');
maba_save_nifti_data(filename_save_nii, results, results_labels, mask_metadata, 'yes');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's save the results in MAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% fprintf('--------------------------------\n');
% fprintf(' SAVING ENCODING AS .MAT        \n');
% fprintf('--------------------------------\n');
% 
% filename_save_mat=strcat(filename_save,'.mat');
% fprintf('Saving MAT file: [\b %s]\b\n', filename_save_mat);
% save(filename_save_mat);

