%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Script by Giacomo Handjaras                                             %%%%%
%%%%% This script create an encoding models from behavioral affective ratings %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
clc

addpath('additional_functions/');
addpath('external_functions/spm12/');

fMRI_run_durations = [ 272 407 402 403 410 390 390 413 ] ; %in TR
fMRI_frequency=0.5; %in Hz
fMRI_TR=1/fMRI_frequency;
%%%there was a fade_in/fade_out in the fMRI scanning session to account for dummy volumes: see the manuscript for more details
fMRI_run_definitions= [ 1 1 2 2 3 3 4 4 ]; %each one of the four sections of the movie during the behavioral experiment was splitted in two fMRI runs

current_folder='../REELMO/jojo-rabbit/emotion-features/subject-avg/';
runs = dir(strcat(current_folder,"*_run-*_emotion-features.h5"));
features=20; %the number of emotion features collected

emotion_model=[];

for r=1:numel(runs)
    INPUT_file=strcat(current_folder, runs(r).name);
    fprintf("    Opening file h5: %s\n",INPUT_file);
    Ratings=h5read(INPUT_file,"/Ratings");
    Categories=h5read(INPUT_file,"/TaggingCategories");
    TempRatings = permute(Ratings,[2 1]);
    SamplingFreq=h5read(INPUT_file,"/RatingsSamplingFrequency");
    DownSamplingFactor=SamplingFreq/fMRI_frequency;
    fMRI_run_mask=find(fMRI_run_definitions==r);

    switch r
        case numel(runs)  %% last run
            while (size(TempRatings,1)<sum(fMRI_run_durations(fMRI_run_mask))*DownSamplingFactor)
                TempRatings=cat(1,TempRatings,zeros(1,features));
            end
        case 1  %% first run
            while (size(TempRatings,1)<sum(fMRI_run_durations(fMRI_run_mask))*DownSamplingFactor)
                TempRatings=cat(1,zeros(1,features),TempRatings);
            end
        otherwise
            while (size(TempRatings,1)<sum(fMRI_run_durations(fMRI_run_mask))*DownSamplingFactor)
                TempRatings=cat(1,TempRatings,zeros(1,features));
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
        emotion_model=cat(1,emotion_model,TempRatingsDS);
        starting_tp=ending_tp+1;
    end

    clear Ratings INPUT_file SamplingFreq DownSamplingFactor starting_tp ending_tp fr fMRI_run_mask TempRatingsDS TempRatings
end


%% plot the model
figure()
imagesc(emotion_model');
yticks([1:numel(Categories)])
yticklabels(Categories);


%% now convolve the model with HRF
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

emotion_model_conv=[];
for i=1:features
emotion_model_conv(:,i)=conv(emotion_model(:,i),hrf);
end

%normalize the data
emotion_model_conv=emotion_model_conv(1:sum(fMRI_run_durations),:);
emotion_model_conv=emotion_model_conv./max(emotion_model_conv(:));

%%Let's save the model in a CSV file
emotion_model_conv_tbl=array2table(emotion_model_conv);
emotion_model_conv_tbl=renamevars(emotion_model_conv_tbl,emotion_model_conv_tbl.Properties.VariableNames, Categories);
writetable(emotion_model_conv_tbl,'models/emotion_model.csv');


