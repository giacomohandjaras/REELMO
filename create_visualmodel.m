%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Script by Giacomo Handjaras                                          %%%%%
%%%%% This script create an encoding models from basic visual features     %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
clc
addpath('additional_functions/');
addpath('external_functions/spm12/');

fMRI_run_durations = [ 272 407 402 403 410 390 390 413 ] ; %in TR
fMRI_frequency=0.5; %in Hz
fMRI_TR=1/fMRI_frequency;
%%%there was a fade_in/fade_out in the fMRI scanning session to account for dummy volumes: see the manuscript for more details
fMRI_run_definitions= [ 1 1 2 2 3 3 4 4 ]; %each one of the four sections of the movie during the behavioral experiment was splitted in two fMRI runs

current_folder='../REELMO/jojo-rabbit/visual-features/';
runs = dir(strcat(current_folder,"*_run-*_global-visual.h5"));
features=4; %the number of visual features collected

model=[];

for r=1:numel(runs)
    INPUT_file=strcat(current_folder, runs(r).name);
    fprintf("    Opening file h5: %s\n",INPUT_file);
    Ratings=[];
    Ratings(:,1)=h5read(INPUT_file,"/contrast");
    Ratings(:,2)=h5read(INPUT_file,"/lightness");
    Ratings(:,3)=h5read(INPUT_file,"/chroma");
    Ratings(:,4)=h5read(INPUT_file,"/hue");
    TempRatings = Ratings; 
    SamplingFreq=h5read(INPUT_file,"/SamplingFrequency");
    DownSamplingFactor=SamplingFreq/fMRI_frequency;
    fMRI_run_mask=find(fMRI_run_definitions==r);

    switch r
        case numel(runs) %% last run
            while (size(TempRatings,1)<sum(fMRI_run_durations(fMRI_run_mask))*DownSamplingFactor)
                TempRatings=cat(1,TempRatings,repmat(min(TempRatings(:)),1,features));
            end
        case 1 %% first run
            while (size(TempRatings,1)<sum(fMRI_run_durations(fMRI_run_mask))*DownSamplingFactor)
                TempRatings=cat(1,repmat(min(TempRatings(:)),1,features),TempRatings);
            end
        otherwise
            while (size(TempRatings,1)<sum(fMRI_run_durations(fMRI_run_mask))*DownSamplingFactor)
                TempRatings=cat(1,TempRatings,repmat(min(TempRatings(:)),1,features));
            end
    end

    starting_tp=1;
    for fr=1:numel(fMRI_run_mask)
        ending_tp=starting_tp-1+(fMRI_run_durations(fMRI_run_mask(fr))*DownSamplingFactor);
        if ending_tp > size(TempRatings,1)
            ending_tp=size(TempRatings,1);
        end
        TempRatingsDS = single(downsample(double( TempRatings(starting_tp:ending_tp,:)),DownSamplingFactor));
        TempRatingsDS(isnan(TempRatingsDS))=0;
        fprintf("    Duration in fMRI: %d\n",size(TempRatingsDS,1));
        model=cat(1,model,TempRatingsDS);
        starting_tp=ending_tp+1;
    end

    clear Ratings INPUT_file SamplingFreq DownSamplingFactor starting_tp ending_tp fr fMRI_run_mask TempRatingsDS TempRatings
end


%% extract features
model_cleaned = filloutliers(model,"linear");
model_cleaned = double(model_cleaned./max(model_cleaned));
visual_energy=model_cleaned;
plot(visual_energy);


%%now convolve the model with HRF
convolve_basis='GAM2'; %%%%% 'GAM' or 'GAM2' for spm12 double gamma
spm_hrf_parameters = [6 16 1 1 6 0 32];
spm_hrf_parameters_t = 16; %%%%% Micro Time

disp(sprintf('HRF function: [\b GAM2]\b'));
hrf = spm_hrf(fMRI_TR,spm_hrf_parameters,spm_hrf_parameters_t);
figure();
plot(hrf); hold on
xlabel('time (in TRs)');
ylabel('relative intensity')
title('GAM2 HRF');
axis square
hold off;
drawnow


visual_energy_conv=[];
for f=1:features
visual_energy_conv(:,f)=conv(visual_energy(:,f),hrf);
end

visual_energy_conv=visual_energy_conv(1:sum(fMRI_run_durations),:);

%%Let's save the model in a CSV file
visual_energy_conv_tbl=array2table(visual_energy_conv);
visual_energy_conv_tbl=renamevars(visual_energy_conv_tbl,visual_energy_conv_tbl.Properties.VariableNames, {'contrast', 'lightness', 'chroma', 'hue'});
writetable(visual_energy_conv_tbl,'models/visual_model.csv');



