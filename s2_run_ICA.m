%% Preprocessing Pipeline Part 2
%
%  -- Bowen Xiu -- 2023

clear all
close all

t_start = tic; % start timer for script
eeglab nogui % starts eeglab

%% variables
% file paths
load_path = char(strcat(pwd,"\eeg\2_preprocessed\"));
save_path = char(strcat(pwd,"\eeg\3_post_ica\"));


% subject file names
all_sub = dir(fullfile(load_path, '*.set'));
all_sub = {all_sub(:).name}';

fprintf("Variables initialized...\n");


%% process files
for n_sub = 1:length(all_sub)
    
    sub_id = all_sub{n_sub}; % gets subject
    sub_start = tic; % start sub-timer
    
    %% load EEG
    % creates EEG file name
    file_name = strcat(load_path, sub_id);  % creates file name for loading
    tmp_idx = find(sub_id == '_');
    sub_id = sub_id(1:(tmp_idx - 1));       % isolates subject id for naming purposes
    set_name = sub_id;

    fprintf("\nCurrent file: %s\n", file_name);

    % load set file
    EEG = pop_loadset('filename', char(file_name));
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);

    %% run ICA
    % calculate number of channels
    n_chan = EEG.nbchan;

    % run ICA
    EEG = pop_runica(EEG, 'extended',1,'interupt','on', 'chanind', 1:n_chan);

    %% save datasets
    % overwrites dataset
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', char(set_name), 'overwrite', 'on', 'gui', 'off');
    % save dataset
    save_name = char(strcat(set_name,"_post_ICA.set"));
    EEG = pop_saveset(EEG, 'filename', save_name, 'filepath', save_path);
    
    %% post-processing steps
    % stop sub-timer
    sub_end = toc(sub_start);
    % print time taken and subject completed
    fprintf("-- Completed: %s\n", sub_id);
    fprintf("-- Elapsed Time: %.2f minutes\n", sub_end/60);

    % clear variables to save memory
    STUDY = [];
    CURRENTSTUDY = 0;
    ALLEEG = [];
    EEG = [];
    CURRENTSET = [];
end

t_end = toc(t_start); % stop timer
fprintf("\n\n\nTotal Time: %.2f minutes\n", t_end/60);
fprintf("--- Goodbye ---\n\n\n");

clear all
close all