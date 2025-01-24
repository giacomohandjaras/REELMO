%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Script by Giacomo Handjaras                                             %%%%%
%%%%% This script create an encoding models from basic acoustic features      %%%%% 
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

current_folder='../REELMO/jojo-rabbit/acoustic-features/';
runs = dir(strcat(current_folder,"*_run-*_power-spectrum.h5"));
features=449; %the number of acoustic spectral features collected

model=[];

for r=1:numel(runs)
    INPUT_file=strcat(current_folder, runs(r).name);
    fprintf("    Opening file h5: %s\n",INPUT_file);
    Ratings=h5read(INPUT_file,"/power");
    Categories=h5read(INPUT_file,"/frequency_bands");
    TempRatings = permute(Ratings,[2 1]);
    SamplingFreq=h5read(INPUT_file,"/SamplingFrequency");
    DownSamplingFactor=SamplingFreq/fMRI_frequency;
    fMRI_run_mask=find(fMRI_run_definitions==r);

    switch r
        case numel(runs)  %% last run
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


%% extract sound energy
sound_energy=sum(model,2);
plot(sound_energy);


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


sound_energy_conv=conv(sound_energy,hrf);

%normalize the data
sound_energy_conv=sound_energy_conv(1:sum(fMRI_run_durations),:);
max_intensity=sound_energy_conv(:);
max_intensity(max_intensity==0)=[];
max_intensity(isnan(max_intensity))=[];
max_intensity(isinf(max_intensity))=[];

%%Let's save the model in a CSV file
sound_energy_conv=double(sound_energy_conv./abs(max(max_intensity)));
sound_energy_conv_tbl=array2table(sound_energy_conv);
sound_energy_conv_tbl=renamevars(sound_energy_conv_tbl,sound_energy_conv_tbl.Properties.VariableNames, {'sound_energy'});
writetable(sound_energy_conv_tbl,'models/acoustic_model.csv');



