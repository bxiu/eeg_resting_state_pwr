%% CALCULATE PSD and Power
%
%  -- Bowen Xiu -- 2023

clear all
close all


fprintf("\n== PSD CALCULATION ==\n");

file_path = strcat(pwd,"\output\eeg_data\");
param_path = strcat(pwd,"\output\subject_parameters\");
filename = uigetfile('*.mat', 'Select a dataset', file_path);
[~, set_name, ~] = fileparts(filename);

data_name = strcat(file_path, filename);
param_name = strcat(param_path, set_name, "_param.mat");

method = input('Select method for calculating the power spectral density (1) pwelch (2) fft: ');

total_st = tic; % start timer for script

fprintf("Loading data...\n");
load(data_name);
fprintf("Data loaded...\n");

fprintf("Loading subject parameters...\n");
load(param_name);
fprintf("Data loaded...\n");

if method == 1
    % calculate PSD using pwelch
    pwelch_method = input('Select pwelch spectrum type (1) ''psd'' (2) ''power'': ');
    
    if pwelch_method == 1
        spectrumtype = 'psd';           % select whether to output the PSD ('psd') or estimate of power at each frequency ('power')
    elseif pwelch_method == 2
        spectrumtype = 'power';
    else
        disp('Type in an existing option.');
    end
    
    % settings
    x = resting_eeg.data';              % input signal
    window = hann(resting_eeg.srate);   % window vector that is 1 sec == EEG srate in length; currenting using Hanning window
    noverlap = resting_eeg.srate/2;     % number of overlapping samples; 50% overlap = EEG srate/2
    psd_res = 0.5;                      % sets the resolution wanted
    nfft = resting_eeg.srate/psd_res;   % number of points in discrete Fourier transform; the 'resolution' == window length/nfft; if nfft is longer than the window
                                        % length, it will be zero-padded
    fs = resting_eeg.srate;             % sampling rate in Hz
    

    [PSD_data, PSD_freqs] = pwelch(x, window, noverlap, nfft, fs, spectrumtype);
    % to convert microvolt^2/Hz to log, do 10*log10(PSD_data)
    % mean_PSD = mean(10*log10(PSD_data), 2);

elseif method == 2
    % calculate PSD using fft
    % settings
    x = resting_eeg.data';              % input signal
    N = length(x);                      % length of signal
    psd_res = 0.5;                      % sets the resolution wanted
    nfft = resting_eeg.srate/psd_res;   % number of points in discrete Fourier transform; the 'resolution' == window length/nfft; if nfft is longer than the window
                                        % length, it will be zero-padded
    fs = resting_eeg.srate;             % sampling rate in Hz  
    
    fft_x = fft(x, nfft, 1).^2;
    PSD_data = (abs(fft_x(1:nfft/2 + 1, :)).^2/N)/fs;
    PSD_freqs = fs/2*linspace(0,1,nfft/2+1);
end

mean_PSD = mean(PSD_data, 2);
iaf_pwr = max(mean_PSD(find(PSD_freqs == 8):find(PSD_freqs == 13)));
iaf = PSD_freqs(mean_PSD == iaf_pwr);

subject_psd = subject;
subject_psd.PSD_data = PSD_data;
subject_psd.mean_PSD = mean_PSD;
subject_psd.iaf = iaf;

fig_psd = figure;
hold on;
plot(PSD_freqs, mean_PSD);
plot(iaf, iaf_pwr,'*','MarkerSize',5);
xlim([subject.filter(1) subject.filter(2)]);
iaf_label = strcat("Mean Individual Alpha Frequency Across Channels: ", num2str(iaf));
title(iaf_label)
if method == 1
    if strcmp(spectrumtype, 'psd')
        ylabel("Power Spectral Density (\muV^2/Hz)");
    elseif strcmp(spectrumtype, 'power')
        ylabel("Power (dB/Hz)");
    end
end
xlabel("Frequency (Hz)");

fig_path = strcat(pwd,"\output\figures\");

if method == 1
    method_name = '_pwelch';
else
    method_name = '_fft';
end

% save figure as .png and .fig
fig_name = strcat(char(fig_path), set_name, method_name, "_psd.png");
saveas(fig_psd, fig_name)
fig_name = strcat(char(fig_path), set_name, method_name, "_psd.fig");
saveas(fig_psd, fig_name)


%% calculate relative and absolute powers
% delta = 1 to 4, theta = 4 to 8, alpha = 8 to 13, beta = 13 to 30
deltaIdx = find(PSD_freqs > 1 & PSD_freqs <= 4);
thetaIdx = find(PSD_freqs > 4 & PSD_freqs <= 8);
alphaIdx = find(PSD_freqs > 8 & PSD_freqs <= 13);
betaIdx  = find(PSD_freqs > 13 & PSD_freqs <= 30);

% compute absolute power for all channels
if method == 1 && strcmp(spectrumtype, 'psd')
    % find integral under the regions of the graph using Simpson's Rule
    % uses simps.m script by Damien Garcia, 2008
    deltaPower = simps(PSD_data(deltaIdx, :));
    thetaPower = simps(PSD_data(thetaIdx, :));
    alphaPower = simps(PSD_data(alphaIdx, :));
    betaPower  = simps(PSD_data(betaIdx, :));
        
elseif method == 2 || strcmp(spectrumtype, 'power')
    % add power together between the frequency intervals
    deltaPower = sum(PSD_data(deltaIdx, :), 1);
    thetaPower = sum(PSD_data(thetaIdx, :), 1);
    alphaPower = sum(PSD_data(alphaIdx, :), 1);
    betaPower  = sum(PSD_data(betaIdx, :), 1);
    
end

% compute relative power
totalPower = deltaPower + thetaPower + alphaPower + betaPower;
delta_relative = deltaPower./totalPower;
theta_relative = thetaPower./totalPower;
alpha_relative = alphaPower./totalPower;
beta_relative  = betaPower./totalPower;

subject_psd.abs_pwr = table(deltaPower', thetaPower', alphaPower', betaPower',...
                            'VariableNames', {'delta', 'theta', 'alpha', 'beta'},...
                            'RowNames', string(subject_psd.ch_names));

subject_psd.rel_pwr = table(delta_relative', theta_relative', alpha_relative', beta_relative',...
                            'VariableNames', {'delta', 'theta', 'alpha', 'beta'},...
                            'RowNames', string(subject_psd.ch_names));

% subject_psd.abs_pwr = struct2table(struct('channel', subject_psd.ch_names,...
%                                           'delta', num2cell(deltaPower),...
%                                           'theta', num2cell(thetaPower),...
%                                           'alpha', num2cell(alphaPower),...
%                                           'beta',  num2cell(betaPower)));
% 
% subject_psd.rel_pwr = struct2table(struct('channel', subject_psd.ch_names,...
%                                           'delta', num2cell(delta_relative),...
%                                           'theta', num2cell(theta_relative),...
%                                           'alpha', num2cell(alpha_relative),...
%                                           'beta',  num2cell(beta_relative)));


abs_fz = subject_psd.abs_pwr{'FZ', :};
abs_cz = subject_psd.abs_pwr{'CZ', :};
abs_pz = subject_psd.abs_pwr{'PZ', :};

rel_fz = subject_psd.rel_pwr{'FZ', :};
rel_cz = subject_psd.rel_pwr{'CZ', :};
rel_pz = subject_psd.rel_pwr{'PZ', :};

output_path = strcat(pwd,"\output\psd\");
output_name = strcat(char(output_path), set_name, "_psd.mat");

save(output_name, 'subject_psd');
fprintf(strcat("Data extracted to \\output\\psd\\ \n"));

% plot absolute power at Fz, Cz, Pz
figure;
subplot(3,1,1);
bar(abs_fz)
set(gca,'xticklabel',{'delta', 'theta', 'alpha', 'beta'});
title('Absolute Power at Fz');
subplot(3,1,2);
bar(abs_cz)
set(gca,'xticklabel',{'delta', 'theta', 'alpha', 'beta'});
title('Absolute Power at Cz');
subplot(3,1,3);
bar(abs_pz)
set(gca,'xticklabel',{'delta', 'theta', 'alpha', 'beta'});
title('Absolute Power at Pz');

% save figure as .png and .fig
fig_name = strcat(char(fig_path), set_name, method_name, "_abs_bandpwr.png");
saveas(fig_psd, fig_name)
fig_name = strcat(char(fig_path), set_name, method_name, "_abs_bandpwr.fig");
saveas(fig_psd, fig_name)

% plot relative power at Fz, Cz, Pz
figure;
subplot(3,1,1);
bar(rel_fz)
set(gca,'xticklabel',{'delta', 'theta', 'alpha', 'beta'});
title('Relative Power at Fz');
subplot(3,1,2);
bar(rel_cz)
set(gca,'xticklabel',{'delta', 'theta', 'alpha', 'beta'});
title('Relative Power at Cz');
subplot(3,1,3);
bar(rel_pz)
set(gca,'xticklabel',{'delta', 'theta', 'alpha', 'beta'});
title('Relative Power at Pz');

% save figure as .png and .fig
fig_name = strcat(char(fig_path), set_name, method_name, "_rel_bandpwr.png");
saveas(fig_psd, fig_name)
fig_name = strcat(char(fig_path), set_name, method_name, "_rel_bandpwr.fig");
saveas(fig_psd, fig_name)

total_end = toc(total_st);
fprintf("\n\n\nTOTAL TIME: %.2f minutes\n", total_end/60);
fprintf("--- Goodbye, friend ---\n\n");
