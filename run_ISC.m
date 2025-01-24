%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Script by Giacomo Handjaras                                          %%%%%
%%%%% This script performs Inter-Subject Correlation (ISC) using           %%%%% 
%%%%% subjects pertaining to the same group.                               %%%%%
%%%%% It requires a mask of grey matter to isolate voxels of interest,     %%%%%
%%%%% additionally a mask of a region of interest to check the quality     %%%%%
%%%%% of each pairing of subjects (within each run of the stimulation),    %%%%%
%%%%% and fMRI data of the subjects.                                       %%%%%
%%%%%                                                                      %%%%%
%%%%% Statistical analyses are performed by using:                         %%%%%
%%%%% 'none' => just reporting the effect size,                            %%%%%
%%%%% 'parametric' => a ttest on r-to-z coeffs,                            %%%%%
%%%%% 'non-parametric' => a wilcoxon rank sum test on r coeffs,            %%%%% 
%%%%% 'permutation' => a non-parametric permutation using a permutation    %%%%%
%%%%%                  schema on fMRI signals, defined by the user.        %%%%%
%%%%%                  P-vals are calculated using Pareto right tail       %%%%%
%%%%%                  approximation.                                      %%%%%
%%%%% Correction for multiple comparisons is handled by using FDR and/or   %%%%%
%%%%% non-parametric permutation using FWE correction.                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

addpath('additional_functions/');
addpath('external_functions/');
addpath('external_functions/NIfTI_tools/');

%%%%% voxels to be analyzed (a binary mask)
filename_mask='../REELMO/jojo-rabbit/fmri/derivatives/group/brain_mask.nii.gz'; 

%%%%% perform a quality check on ISC in a region of interest (potentially with a decent SNR)
quality_check='no'; 
filename_mask_qc='../REELMO/jojo-rabbit/fmri/derivatives/group/V1.nii.gz'; %%% Mask of interest to perform quality check of ISC across subjects

%%%%% files to be saved
filename_save='../REELMO/jojo-rabbit/fmri/derivatives/group/Results_ISC_parametric'; %%% a nifti .nii and a matlab .mat will be created at the end of the script

%%%%% fMRI subjects ID
groupA={ '01','02','03','04','05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20' };
groupA_path='../REELMO/jojo-rabbit/fmri/derivatives/'; %%% the absolute path (in BIDS format)
groupA_file='_epi_cleaned_reml_2mni.nii.gz'; %%% the actual fMRI preprocessed files to load (in BIDS format) 

%%%%% the duration of each run in volumes. Actually the run_durations is summed to obtain the total duration of the timeserie
run_durations=[ 272 407 402 403 410 390 390 413 ] ; %%%run duration in TRs

%%%%% censored time points
run_censored_timepoints=[]; %%% if someone wants to remove specific timepoins

%%%%% perform statistics! It could be 'none' if we aim to obtain only the effect size, 'parametric' if we want to use a ttest on r-to-z coeffs, 'non-parametric' using a wilcoxon rank sum test, 'permutation' using a non-parametric permutation test 
statistics='parametric';
statistics_correction='fdr95'; %%%%%'fdr95' or 'fdr01'
statistics_permutation_data=''; %%% permutation schema for statistic='permutation'

%%%%% define the precision of the data
data_precision='single';

%%%%% number of CPU cores to be used
CPUs=6;

%%%%% number of voxel bins (to limit memory use and speedup computations)
voxels_bins=25;

%%%%% number in % below which a voxel is discarded from analyses (e.g., not all the voxels are acquired in all participants)
tollerance=50; %%%% 50, means that the voxel is included in the analysis only when at least 50% of the participants have data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's prepare the environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------\n');
fprintf(' PREPARING ENVIRONMENT          \n');
fprintf('--------------------------------\n');

%%%%% set the flag for the quality check
quality_check_flag = sum(ismember(lower(quality_check),{'yes','y'}))>0;

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
statistics_flag = sum(ismember(statistics,{'none','parametric','non-parametric','permutation'}))>0;
if statistics_flag
    fprintf('Statistics defined: [\b %s]\b\n', statistics);
    groupA_random_schema=[];
    groupA_permutations=[];
    groupA_permuted_subjects=[];
else
    error('Statistics "[\b %s]\b" was not supported!', statistics);
end

%%%%% open permutation schema
if strcmp(statistics,'permutation')
    fprintf('Opening permutation schema: [\b %s]\b\n', statistics_permutation_data);
    groupA_random_schema=load(statistics_permutation_data);
    groupA_permutations=size(groupA_random_schema.permutation_tps,1);
    fprintf('Number of permutations: %d\n', groupA_permutations);
    groupA_permuted_subjects=size(groupA_random_schema.permutation_tps,2);
    fprintf('Number of permuted subjects supported: %d\n', groupA_permuted_subjects);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's open the fMRI masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------\n');
fprintf(' LOADING MASKS                  \n');
fprintf('--------------------------------\n');

[mask_volume,mask_metadata]=maba_load_nifti_mask(filename_mask,'yes');
if quality_check_flag
    [mask_qc_volume,mask_qc_metadata]=maba_load_nifti_mask(filename_mask_qc,'yes');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's open the subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------\n');
fprintf(' LOADING fMRI                   \n');
fprintf('--------------------------------\n');

groupA_subjects=numel(groupA);
groupA_pairings=groupA_subjects*(groupA_subjects-1)/2;
groupA_tps=sum(run_durations);
groupA_alldata=nan(groupA_subjects,mask_metadata.voxels,groupA_tps,data_precision);
if quality_check_flag
    groupA_qc_alldata=nan(groupA_subjects,mask_qc_metadata.voxels,groupA_tps,data_precision);
end

for s=1:groupA_subjects
    groupA_filename=strcat(groupA_path,'sub-',groupA{s},'/','sub-',groupA{s},groupA_file);
    [groupA_volume]=maba_load_nifti_data(groupA_filename,'no');
    groupA_data=maba_get_nifti_data(groupA_volume,mask_metadata.coordinates);
    groupA_alldata(s,:,:)=groupA_data;
    if quality_check_flag
        groupA_qc_data=maba_get_nifti_data(groupA_volume,mask_qc_metadata.coordinates);
        groupA_qc_alldata(s,:,:)=groupA_qc_data;
        clear  groupA_qc_data
    end
    clear groupA_data  groupA_volume groupA_filename
end


%%%%% check precision again
if strcmp(data_precision,'single')
    groupA_alldata=single(groupA_alldata);
    if quality_check_flag
        groupA_qc_alldata=single(groupA_qc_alldata);
    end
else
    groupA_alldata=double(groupA_alldata);
    if quality_check_flag
        groupA_qc_alldata=double(groupA_qc_alldata);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Before moving on, let's check the goodness of data by exploring ISC within filename_mask_qc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if quality_check_flag

    fprintf('--------------------------------\n');
    fprintf(' QUALITY CHECK                  \n');
    fprintf('--------------------------------\n');

    groupA_qc_alldata_mad=mad(groupA_qc_alldata,0,3); %%% let's extract MAD. If it is zero means no signal (mask outside fmri data)!
    groupA_qc_alldata_mask=groupA_qc_alldata_mad>0;
    groupA_qc_alldata_avg=nan(groupA_subjects,groupA_tps,data_precision);

    for s=1:groupA_subjects
        groupA_qc_alldata_avg(s,:)=nanmean(groupA_qc_alldata(s,groupA_qc_alldata_mask(s,:),:),2);
    end

    %%%%% censor data
    groupA_qc_alldata_avg(:,run_censored_timepoints)=[];
    groupA_qc_ISC=1-pdist(groupA_qc_alldata_avg,'correlation');
    groupA_qc_ISC_sq=squareform(groupA_qc_ISC);
    %groupA_qc_ISC_sq(groupA_qc_ISC_sq<prctile(groupA_qc_ISC,5))=nan;
    groupA_qc_ISC_sq(groupA_qc_ISC_sq<=0)=nan;
    groupA_qc_ISC_tree = linkage(1-groupA_qc_ISC,'average'); 
    groupA_qc_ISC_order = optimalleaforder(groupA_qc_ISC_tree,1-groupA_qc_ISC);

    groupA_qc_ISC_run=nan(numel(run_durations),numel(groupA_qc_ISC));
    for r=1:numel(run_durations)
        groupA_qc_alldata_avg_run=groupA_qc_alldata_avg(:,run_tps==r);
        groupA_qc_ISC_run(r,:)=1-pdist(groupA_qc_alldata_avg_run,'correlation');
    end
    groupA_qc_ISC_run_mask=groupA_qc_ISC_run<prctile(groupA_qc_ISC_run(:),0.25);
    groupA_qc_ISC_run_mask=sum(groupA_qc_ISC_run_mask)>0;
    clear r groupA_qc_alldata_avg_run

    labels_isc_temp=cell(numel(groupA),numel(groupA));
    for r=1:numel(groupA)
    for c=1:numel(groupA)
        labels_isc_temp(r,c)=strcat(groupA(r),'~',groupA(c));
    end
    end
    labels_isc_indx=squareform(1:numel(groupA_qc_ISC));
    labels_isc=cell(numel(groupA_qc_ISC),1);
    for r=1:numel(groupA_qc_ISC)
        labels_isc_temp_temp=labels_isc_temp(labels_isc_indx==r);
        labels_isc(r)=labels_isc_temp_temp(1);
    end
    clear labels_isc_temp labels_isc_temp_temp labels_isc_indx r c


    %%%%% let's draw three figures: 
    % a1) an histogram with ISC values, 
    % a2) squareform of ISC (done using quick and dirty pdist)
    % b) squareform of ISC, after sorting of the subjects according to their similarities
    % c) check ISC in each run individually
    figure();
    subplot(1,2,1);
        histogram(groupA_qc_ISC);
        xlabel('ISC (corr coeff)');
        ylabel('subject''s pairings');
        axis square;
        title(sprintf('Average ISC: %.2f',nanmean(groupA_qc_ISC)));
    subplot(1,2,2);
        imagesc(groupA_qc_ISC_sq);
        %imagesc(groupA_qc_ISC_sq(groupA_qc_ISC_order,groupA_qc_ISC_order));
        caxis([prctile(groupA_qc_ISC,5),prctile(groupA_qc_ISC,95)]);
        quality_check_colormap=colormap('Autumn');
        quality_check_colormap(1,:)=[0,0,0];
        colormap(quality_check_colormap);
        colorbar;
        xlabel('Subjects');
        ylabel('Subjects');
        xticks(1:groupA_subjects);
        yticks(1:groupA_subjects);
        xticklabels(groupA);
        yticklabels(groupA);
        %xticklabels(groupA(groupA_qc_ISC_order));
        %yticklabels(groupA(groupA_qc_ISC_order));
        xtickangle(45);
        axis square;
        title(sprintf('Individual ISC'));
    sgtitle('Quality Check of ISC');

    figure();
    imagesc(groupA_qc_ISC_sq(groupA_qc_ISC_order,groupA_qc_ISC_order));
    caxis([prctile(groupA_qc_ISC,5),prctile(groupA_qc_ISC,95)]);
    quality_check_colormap=colormap('Autumn');
    quality_check_colormap(1,:)=[0,0,0];
    colormap(quality_check_colormap);
    colorbar;
    xlabel('Subjects');
    ylabel('Subjects');
    xticks(1:groupA_subjects);
    yticks(1:groupA_subjects);
    xticklabels(groupA(groupA_qc_ISC_order));
    yticklabels(groupA(groupA_qc_ISC_order));
    xtickangle(45);
    axis square;
    title(sprintf('Individual ISC, reordered'));

    figure();
    %%%%% plot the average effect across runs and a linear fit
    scatter(1:numel(run_durations),mean(groupA_qc_ISC_run,2),60,[0,0,0],'filled','s'); hold on
    ls=lsline;
    ls.LineWidth=3;
    ls.LineStyle='--';
    ls.Color=[0,0,0];

    %%%%% plot the zero
    plot(0:numel(run_durations)+1,repmat([0],numel(run_durations)+2,1),'Color','k','LineStyle','--','LineWidth',0.5); hold on

    %%%%% plot the pairings, changing colors as a function of their values
    for s=1:numel(groupA_qc_ISC)
        if(groupA_qc_ISC_run_mask(s)==1)
        random_color=rand(1,3);
        random_alpha=1;
        random_size=20;
        random_jitter=0;
        else
        random_color=[0.5,0.5,0.5];
        random_alpha=0.3;
        random_size=10;
        random_jitter=randn(1,numel(run_durations))/8;
        end    
        scatter([1:numel(run_durations)]+random_jitter,groupA_qc_ISC_run(:,s),random_size,random_color,'filled','MarkerFaceAlpha',random_alpha);

        if(groupA_qc_ISC_run_mask(s)==1)
        plot([1:numel(run_durations)]+random_jitter,groupA_qc_ISC_run(:,s),'Color',random_color,'LineWidth',2);
        end
    end
    %%%%% write the text (at the end to avoid overalp)    
    for s=1:numel(groupA_qc_ISC)
        if(groupA_qc_ISC_run_mask(s)==1)
        text([1:numel(run_durations)]+random_jitter,groupA_qc_ISC_run(:,s),repmat(labels_isc(s),numel(run_durations),1),'FontSize',6);
        end 
    end
    axis square;
    pbaspect([3 1 1])
    xlabel('Run');
    xticks([1:numel(run_durations)])
    ylabel('correlation coeff');
    title('ISC across pairings and runs');
    hold off
    drawnow
    input('Press ''Enter'' to start calculation of ISC...','s');
    clear s random_color random_alpha random_size random_jitter ls
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's measure ISC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------\n');
fprintf(' EVALUATE ISC                   \n');
fprintf('--------------------------------\n');

%%%%% To reduce memory usage with parallelization (parfor), we need to split
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
clear starting_point ending_point

%%%%% Now perform a for loop for each bin
if CPUs>1
    CPUs_handle=parpool('local',CPUs);
end

groupA_results=nan(mask_metadata.voxels,groupA_pairings);
groupA_results_pvals=nan(mask_metadata.voxels,1);

for b=1:voxels_bins
    tic
    fprintf('bin %d out of %d....\n',b,voxels_bins);
    current_groupA_data=groupA_alldata(:,mask_voxels_bins(b,1):mask_voxels_bins(b,2),:);
    [current_groupA_results,current_groupA_results_pvals,~]=maba_measure_ISC_onegroup(current_groupA_data,'yes',tollerance,run_censored_timepoints,statistics,groupA_random_schema);
    groupA_results(mask_voxels_bins(b,1):mask_voxels_bins(b,2),:)=current_groupA_results;
    groupA_results_pvals(mask_voxels_bins(b,1):mask_voxels_bins(b,2),:)=current_groupA_results_pvals;
    toc
end

if CPUs>1
    delete(CPUs_handle);
    clear CPUs_handle
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
    results=nanmean(groupA_results,2);
    results_labels={'ISC'};
end

if strcmp(statistics,'parametric')
    groupA_results_effect=nanmean(groupA_results,2);
    groupA_results_pvals_log=abs(log10(groupA_results_pvals));
    groupA_results_pvals_log(isnan(groupA_results_pvals_log))=0;

    groupA_results_pvals_corr=maba_correct_pvals(groupA_results_pvals,'yes',statistics_correction);
    groupA_results_pvals_corr_log=abs(log10(groupA_results_pvals_corr));
    groupA_results_pvals_corr_log(isnan(groupA_results_pvals_corr_log))=0;

    results=cat(2,groupA_results_effect,groupA_results_pvals_log,groupA_results_pvals_corr_log);
    results_labels=strcat('ISC',{' '},['log_pval-raw_',statistics],{' '}, ['log_pval-',statistics_correction,'_', statistics]);
    clear groupA_results_effect groupA_results_pvals_log groupA_results_pvals_corr_log
end

if strcmp(statistics,'non-parametric')
    groupA_results_effect=nanmean(groupA_results,2);
    groupA_results_pvals_log=abs(log10(groupA_results_pvals));
    groupA_results_pvals_log(isnan(groupA_results_pvals_log))=0;

    groupA_results_pvals_corr=maba_correct_pvals(groupA_results_pvals,'yes',statistics_correction);
    groupA_results_pvals_corr_log=abs(log10(groupA_results_pvals_corr));
    groupA_results_pvals_corr_log(isnan(groupA_results_pvals_corr_log))=0;

    results=cat(2,groupA_results_effect,groupA_results_pvals_log,groupA_results_pvals_corr_log);
    results_labels=strcat('ISC',{' '},['log_pval-raw_',statistics],{' '}, ['log_pval-',statistics_correction,'_', statistics]);
    clear groupA_results_effect groupA_results_pvals_log groupA_results_pvals_corr_log
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's save the results in nifti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--------------------------------\n');
fprintf(' SAVING DATA IN 3D          \n');
fprintf('--------------------------------\n');

filename_save_nii=strcat(filename_save,'.nii');
maba_save_nifti_data(filename_save_nii, results, results_labels, mask_metadata, 'yes');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's save the results in MAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('--------------------------------\n');
% fprintf(' SAVING DATA AS .MAT        \n');
% fprintf('--------------------------------\n');
% 
% filename_save_mat=strcat(filename_save,'.mat');
% fprintf('Saving MAT file: [\b %s]\b\n', filename_save_mat);
% save(filename_save_mat);

