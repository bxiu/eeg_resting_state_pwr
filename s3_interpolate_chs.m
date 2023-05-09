%% Preprocessing Pipeline Part 3
%% REMOVE ICA Artifacts before running script
%
%  -- Bowen Xiu -- 2023

% interpolates removed channels
clear all
close all

t_start = tic; % start timer for script
eeglab nogui % starts eeglab

%% variables
% file paths
load_path = char(strcat(pwd,"\eeg\4_artifacts_rmd\"));
save_path = char(strcat(pwd,"\eeg\5_cleaned\"));
save_path_mat = char(strcat(pwd,"\output\subject_parameters\"));
save_path_mat_eeg = char(strcat(pwd,"\output\eeg_data\"));

% subject file names
all_sub = dir(fullfile(load_path, '*.set'));
all_sub = {all_sub(:).name}';

% dataset with all channels
interp_data = 'interpolation_channels.set';

%% load dataset with all channels
EEG = pop_loadset('filename',interp_data,'filepath',pwd);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);

%% process files for selected condition
for n_sub = 1:length(all_sub)
    
    sub_id = all_sub{n_sub};

    %% load EEG
    % creates EEG file name
    file_name = strcat(load_path, sub_id);  % creates file name for loading
    tmp_idx = find(sub_id == '_');
    sub_id = sub_id(1:(tmp_idx - 1));       % isolates subject id for naming purposes
    set_name = sub_id;
    
    param_name = strcat(save_path_mat, set_name, "_param.mat");
    load(param_name);

    fprintf("\nCurrent file: %s\n", file_name);

    % load set file
    EEG = pop_loadset('filename', char(file_name));
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);

    %% interpolate removed channels
    EEG = pop_interp(EEG, ALLEEG(1).chanlocs, 'spherical');
    
    %% average reference
    EEG = pop_reref(EEG, []);
    EEG = eeg_checkset(EEG);

    %% save dataset
    % overwrites dataset
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'setname', char(set_name), 'overwrite', 'on', 'gui', 'off');
    % save dataset
    save_name = char(strcat(set_name,"_cleaned.set"));
    EEG = pop_saveset(EEG, 'filename', save_name, 'filepath', save_path);
    
    % save dataset as .mat file for easier analysis
    data_name = strcat(save_path_mat_eeg, set_name, ".mat");
    resting_eeg.data = EEG.data;
    resting_eeg.srate = EEG.srate;
    save(data_name, 'resting_eeg');
    
    subject.rmd_components = size(EEG.icaweights, 2) - size(EEG.icaweights, 1);
    subject.ch_names = {EEG.chanlocs.labels};
    save(param_name, 'subject');

    % clear variables to save memory
    ALLEEG = pop_delset(ALLEEG, [2]);

end

t_end = toc(t_start); % stop timer
fprintf("\n\n\nTotal Time: %.2f minutes\n", t_end/60);
fprintf("--- Goodbye ---\n\n\n");

clear all
close all