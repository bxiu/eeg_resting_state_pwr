%% Preprocessing Pipeline Part 1
%% Remove extremely bad sections of the EEG data BEFORE RUNNING SCRIPT
%
%  -- Bowen Xiu -- 2023


clear all
close all

fprintf("\n== STEP 2 - PREPROCESS DATA ==\n");

total_st = tic; % start timer for script
eeglab nogui % starts eeglab

%% variables
% file paths
load_path = char(strcat(pwd,"\eeg\1_no_spikes\"));
save_path = char(strcat(pwd,"\eeg\2_preprocessed\"));
save_path_mat = char(strcat(pwd,"\output\subject_parameters\"));

% EEG filter
eeg_low = 0.5;
eeg_high = 40;

% study sampling rate
study_fs = 250;

% subject file names
all_sub = dir(fullfile(load_path, '*.set'));
all_sub = {all_sub(:).name}';

fprintf("Variables initialized...\n");

%% filter each dataset
for n_sub = 1:length(all_sub)
    
    % gets subject id
    sub_id = all_sub{n_sub}; 
    file_name = strcat(load_path, sub_id);  % creates file name for loading
    tmp_idx = find(sub_id == '_');
    sub_id = sub_id(1:(tmp_idx - 1));       % isolates subject id for naming purposes
    
    % create subject parameters file
    subject = struct();
    subject.id = sub_id;
    subject.filter = [eeg_low eeg_high];
    subject.fs = study_fs;

    % process files for selected condition
    t_start = tic;

    % creates EEG dataset name
    set_name = strcat(sub_id, "_preprocessed");

    fprintf("\n\nCurrent file: %s\n", set_name)

    % load file
	% EEG = pop_loadcnt(char(file_name), 'dataformat', 'auto', 'memmapfile', '');
    EEG = pop_loadset('filename',char(file_name));
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', char(set_name), 'gui', 'off');
    
    % Remove channels and add channel locations
    bad_ch = input('Type in each bad channel as chars in a {} bracket. Separate the bad channels with a comma (e.g., {''Cz'', ''C1'', ''Pz''}): ');
    usual_rm_chs = {'CB1' 'CB2' 'VEO' 'HEO' 'EKG' 'EMG'};
    EEG = pop_select(EEG,'nochannel',unique([usual_rm_chs bad_ch]));
    EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp');
    
    subject.bad_chs = bad_ch;
    subject.other_rmd_chs = usual_rm_chs;
    subject.n_bad_chs = length(bad_ch);
    subject.n_rmd_chs = length(unique([usual_rm_chs bad_ch]));

    % filter EEG
    EEG = pop_eegfiltnew(EEG, eeg_low, eeg_high);

    % resamples EEG
    EEG = pop_resample(EEG, study_fs);

    %% Save Dataset
    % overwrites dataset
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', char(set_name), 'overwrite', 'on', 'gui', 'off');
    % save dataset
    save_name = char(strcat(set_name,".set"));
    EEG = pop_saveset(EEG, 'filename', save_name, 'filepath', save_path);
    
    param_name = strcat(save_path_mat, sub_id, "_param.mat");
    save(param_name, 'subject');

    % clear variables to save memory
    STUDY = [];
    CURRENTSTUDY = 0;
    ALLEEG = [];
    EEG = [];
    CURRENTSET = [];

    t_end = toc(t_start);

    fprintf("\n\nElapsed Time: %.2f minutes\n", t_end/60);
end

total_end = toc(total_st);
fprintf("\n\n\nTotal Time: %.2f minutes\n", total_end/60);
fprintf("--- Goodbye ---\n\n");
