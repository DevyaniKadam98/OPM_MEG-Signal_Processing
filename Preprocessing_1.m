% Load OPM-MEG data
data = load('Path to your MEG data....002\OPM\Motor\run02.meg.mat');

% Extract MEG signals (magnetic field values) and convert to nanotesla
meg_signals = data.bexp * 1e9;  % Convert MEG signals from T to nT

% Extract EMG signals and channel names
emg_signals = data.bexp_ext;  % EMG signals in Volts
emg_channel_names = {'trigger_emg', 'sync_pulse', 'trigger_pulse'};

% Extract the trigger pulse channel
trigger_pulse = emg_signals(3, :); % Third EMG channel is 'trigger_pulse'

% Define channel names for MEG data
meg_channel_names = {'M001_y', 'M001_z', 'M002_y', 'M002_z', 'M003_y', 'M003_z', ...
                     'M004_y', 'M004_z', 'M005_y', 'M005_z', 'M006_y', 'M006_z', ...
                     'M007_y', 'M007_z', 'M008_y', 'M008_z', 'M009_y', 'M009_z', ...
                     'M010_y', 'M010_z', 'M011_y', 'M011_z', 'M012_y', 'M012_z', ...
                     'M013_y', 'M013_z', 'M014_y', 'M014_z', 'M015_y', 'M015_z'};

% Define sampling frequency
Fs = 2000; % Sampling frequency in Hz

% Time vector for full duration
time_vector = (0:size(meg_signals, 2)-1) / Fs; % Time vector in seconds

%% Step 1: Detrend the raw MEG signals
meg_signals_detrended = detrend(meg_signals')';  % Remove linear trends from each channel

%% Step 2: Filter the detrended MEG signals
% High-pass filter at 1 Hz
[b_hp, a_hp] = butter(4, 1 / (Fs / 2), 'high');
meg_filtered_hp = filtfilt(b_hp, a_hp, meg_signals_detrended')';

% Low-pass filter at 200 Hz
[b_lp, a_lp] = butter(4, 200 / (Fs / 2), 'low');
meg_filtered_lp = filtfilt(b_lp, a_lp, meg_filtered_hp')';

% Band-stop filter at 60 Hz
[b_bs, a_bs] = butter(4, [56 64] / (Fs / 2), 'stop');
meg_signals_filtered1 = filtfilt(b_bs, a_bs, meg_filtered_lp')';

% Band-stop filter at 120 Hz
[b_bs, a_bs] = butter(4, [116 124] / (Fs / 2), 'stop');
meg_signals_filtered = filtfilt(b_bs, a_bs, meg_signals_filtered1')';

%% Step 3: Apply CAR to the filtered MEG signals
% Calculate the mean across all channels (common average)
mean_signal = mean(meg_signals_filtered, 1);

% Subtract the common average from each channel
meg_signals_car = meg_signals_filtered - mean_signal;

%% Step 4: Plot the Data
figure('Name', 'MEG Signal Processing', 'Color', 'w', 'Position', [0, 0, 1200, 800]);

% Subplot 1: Raw MEG signals
subplot(4, 1, 1);
hold on;
for ch = 1:size(meg_signals, 1)
    plot(time_vector, meg_signals(ch, :), 'DisplayName', meg_channel_names{ch});
end
xlabel('Time (seconds)');
ylabel('Amplitude (nT)');
xlim([0 316]);
title('Raw MEG Signals');
grid on;
hold off;

% Subplot 2: Detrended MEG signals
subplot(4, 1, 2);
hold on;
for ch = 1:size(meg_signals_detrended, 1)
    plot(time_vector, meg_signals_detrended(ch, :), 'DisplayName', meg_channel_names{ch});
end
xlabel('Time (seconds)');
ylabel('Amplitude (nT)');
xlim([0 316]);
title('Detrended MEG Signals');
grid on;
hold off;

% Subplot 3: Filtered MEG signals
subplot(4, 1, 3);
hold on;
for ch = 1:size(meg_signals_filtered, 1)
    plot(time_vector, meg_signals_filtered(ch, :), 'DisplayName', meg_channel_names{ch});
end
xlabel('Time (seconds)');
ylabel('Amplitude (nT)');
xlim([0 316]);
title('Filtered MEG Signals (1 Hz High-pass, 200 Hz Low-pass, 60 Hz Band-stop)');
grid on;
hold off;

% Subplot 4: CAR-corrected MEG signals
subplot(4, 1, 4);
hold on;
for ch = 1:size(meg_signals_car, 1)
    plot(time_vector, meg_signals_car(ch, :), 'DisplayName', meg_channel_names{ch});
end
xlabel('Time (seconds)');
ylabel('Amplitude (nT)');
xlim([0 316]);
title('CAR-Corrected MEG Signals');
grid on;
hold off;

sgtitle('MEG Signal Processing Steps'); % Overall title for the figure

%% Step 5: Save the CAR-corrected data
save('CAR_Corrected_MEG_with_Processed_EMG_002_2.mat', 'meg_signals_car', 'meg_channel_names', 'emg_channel_names','trigger_pulse', 'Fs');
disp('CAR-corrected MEG data saved to "CAR_Corrected_MEG_with_Processed_EMG_002_.mat".');


%% Step 6: Plot Power Spectrum of Raw and Filtered MEG Data
% Parameters for power spectrum
nfft = 2^nextpow2(size(meg_signals, 2)); % FFT length (next power of 2)
freq = (0:nfft/2-1) * Fs / nfft;        % Frequency axis (up to Nyquist frequency)

% Initialize matrices to store power spectra
power_spectrum_raw = zeros(size(meg_signals, 1), nfft/2); % [channels x frequencies]
power_spectrum_filtered = zeros(size(meg_signals, 1), nfft/2); % [channels x frequencies]

% Compute power spectrum for raw and filtered signals for each channel
for ch = 1:size(meg_signals, 1)
    % Raw MEG signals
    fft_raw = fft(meg_signals(ch, :), nfft); % FFT for raw signals
    power_spectrum_raw(ch, :) = abs(fft_raw(1:nfft/2)).^2 / nfft; % Power spectrum (normalize by nfft)

    % Filtered MEG signals
    fft_filtered = fft(meg_signals_filtered(ch, :), nfft); % FFT for filtered signals
    power_spectrum_filtered(ch, :) = abs(fft_filtered(1:nfft/2)).^2 / nfft; % Power spectrum (normalize by nfft)
end

% Plot Power Spectrum for Raw and Filtered Signals
figure('Name', 'Power Spectrum Comparison', 'Color', 'w', 'Position', [0, 0, 1200, 800]);

% Subplot 1: Power spectrum of raw MEG signals
subplot(2, 1, 1);
hold on;
for ch = 1:size(meg_signals, 1)
    plot(freq, 10*log10(power_spectrum_raw(ch, :)), 'DisplayName', meg_channel_names{ch}); % Plot in dB
end
xlabel('Frequency (Hz)');
xlim([0 200]);
ylabel('Power (dB)');
title('Power Spectrum of Raw MEG Signals');
grid on;
legend('show', 'Location', 'northeastoutside'); % Legend outside the plot
hold off;

% Subplot 2: Power spectrum of filtered MEG signals
subplot(2, 1, 2);
hold on;
for ch = 1:size(meg_signals_filtered, 1)
    plot(freq, 10*log10(power_spectrum_filtered(ch, :)), 'DisplayName', meg_channel_names{ch}); % Plot in dB
end
xlabel('Frequency (Hz)');
xlim([0 200]);
ylabel('Power (dB)');
title('Power Spectrum of Filtered MEG Signals');
grid on;
legend('show', 'Location', 'northeastoutside'); % Legend outside the plot
hold off;

