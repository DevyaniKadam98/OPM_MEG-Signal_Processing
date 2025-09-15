% Load the averaged trials data
    data = load('Cut_Trials_MEG_combined_093.mat'); 

    % Extract relevant data
    trials = data.combined_trials; % Trials data: [sensors x time points x trials]
    trial_time_axis = data.combined_trial_time_axis; % Time axis for trials in ms
    meg_channel_names = data.meg_channel_names; % Sensor names
    sf = data.Fs; % Sampling frequency (Hz)

    % Step 1: Average the trials for each sensor
    avg_trials = mean(trials, 3); % Average across the 3rd dimension (trials)

    % Step 2: Baseline correction (remove mean of pre-trigger period)
    pre_trigger_samples = find(trial_time_axis < 0); % Samples before the trigger
    baseline_mean = mean(avg_trials(:, pre_trigger_samples), 2); % Compute baseline mean
    avg_trials_bc = avg_trials - baseline_mean; % Subtract baseline

    % Step 3:Time-Frequency Analysis
    TW = 200; % Window length (samples)
    TS = 100; % Step size (samples)
    plotTimeFrequency(avg_trials_bc, sf, TW, TS, meg_channel_names, 'Time-Frequency Analysis (Averaged Trials)');


function plotTimeFrequency(data, fs, TW, TS, channel_names, title_str, ~)
    % Computes and plots the time-frequency representation of multi-channel data
    % with trigger location marked.
    %
    % Inputs:
    %   data - 2D array (channels x time points)
    %   fs - Sampling frequency (e.g., 2000 Hz)
    %   TW - Window length (in samples)
    %   TS - Step size (in samples)
    %   channel_names - Cell array of channel names
    %   title_str - Title for the overall figure
    %   trial_time_axis - Time axis for trials (in ms)

    % Parameters
    num_channels = size(data, 1); % Number of channels
    num_timepoints = size(data, 2); % Number of time points
    nfft = TW; % Number of FFT points (same as window length)

    % Define the Blackman-Harris window
    winfun = blackmanharris(TW);

    % Define time window start indices
    TWs = TW:TS:num_timepoints;

    % Generate frequency vector in Hz
    freq = (0:(nfft/2)-1) * (fs / nfft);

    % Limit frequency to 100 Hz
    freq_limit = 100; % Limit in Hz
    freq_idx = freq <= freq_limit; % Frequency indices <= 100 Hz
    freq = freq(freq_idx); % Subset frequency vector

    % Initialize power spectrum array
    ps_on = zeros(numel(TWs), numel(freq), num_channels);

    % Iterate over channels for time-frequency analysis
    for ch = 1:num_channels
        disp(['Processing Channel ' num2str(ch) ' / ' num2str(num_channels)]);
        for tstep = 1:numel(TWs)
            % Extract data segment for the current time window
            segment_start = TWs(tstep) - TW + 1;
            segment_end = TWs(tstep);
            if segment_start > 0 && segment_end <= num_timepoints
                dat = data(ch, segment_start:segment_end);
                dat = dat .* winfun'; % Apply Blackman-Harris window
                f = abs(fft(dat, nfft)); % Compute FFT magnitude
                ps_on(tstep, :, ch) = f(freq_idx); % Store frequencies within limit
            end
        end
    end

    % Compute baseline correction
    baseline_ps = mean(ps_on(1:3, :, :), 1); % Mean power in the first 3 time windows
    baseline_corrected_ps = ps_on ./ baseline_ps; % Always divide for freq domain baseline correction

    % Log-transform baseline-corrected power spectrum
    log_baseline_corrected_ps = log10(max(baseline_corrected_ps, eps)); % Add small constant to avoid log(0)

    % Generate time vector in ms
    time_freq_axis = (TWs - TW / 2) / fs * 1000; % Convert to ms

    % Dynamically calculate grid layout
    num_rows = ceil(sqrt(num_channels));
    num_cols = ceil(num_channels / num_rows);
    save('Power_time_freq_values.mat', 'time_freq_axis', 'freq', 'log_baseline_corrected_ps');

    % Plot Time-Frequency Representations
    figure('Name', 'Time-Frequency Plots', 'Color', 'w', 'Position', [0, 0, 1200, 800]);
    for ch = 1:num_channels
        subplot(num_rows, num_cols, ch);

        % Plot the spectrogram
        imagesc(time_freq_axis, freq, squeeze(log_baseline_corrected_ps(:, :, ch))');
        axis tight;
        set(gca, 'YDir', 'normal'); % Correct frequency orientation
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
        title(channel_names{ch});
        colorbar;
        clim([-2 2]); % Adjust color scale

        % Plot the trigger location
        trigger_time = 501; % trigger is at 501 ms 
        hold on;
        plot([trigger_time trigger_time], [min(freq) max(freq)], 'r', 'LineWidth', 1); % Red line for trigger
        hold off;
    end

    sgtitle(title_str); % Add overall title
end
