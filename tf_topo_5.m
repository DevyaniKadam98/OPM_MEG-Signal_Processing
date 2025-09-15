    % Load the averaged trials data
    data = load('Cut_Trials_MEG_combined_002.mat'); 

    % Extract relevant data
    trials = data.combined_trials;           % Trials data: [sensors x time points x trials]
    trial_time_axis = data.combined_trial_time_axis; % Time axis for trials in ms
    meg_channel_names = data.meg_channel_names;      % Sensor names (e.g., 'M001_y', 'M001_z', …)
    sf = data.Fs;                           % Sampling frequency (Hz)

    % Separate Y and Z channels based on naming convention
    y_idx = contains(meg_channel_names, 'y');
    z_idx = contains(meg_channel_names, 'z');

    % Step 1: Average the trials for each sensor
    avg_trials = mean(trials, 3);  % [sensors x time points]
    avg_trials_y = avg_trials(y_idx, :);
    avg_trials_z = avg_trials(z_idx, :);

    % Step 2: Baseline correction (remove mean of pre-trigger period)
    pre_trigger_samples = find(trial_time_axis < 0); % Samples before the trigger
    baseline_mean_y = mean(avg_trials_y(:, pre_trigger_samples), 2);
    baseline_mean_z = mean(avg_trials_z(:, pre_trigger_samples), 2);
    avg_trials_y_bc = avg_trials_y - baseline_mean_y;
    avg_trials_z_bc = avg_trials_z - baseline_mean_z;

    % Step 3: Time-Frequency Analysis parameters
    TW = 200; % Window length (samples)
    TS = 100; % Step size (samples)

    % Define sensor layout (4 rows x 5 columns)
    % Use -1 for empty slots.
    % The layout is arranged with “Nose” at the top.
    %
    %       Nose
    %  [-1, 14, 13, 12, -1]
    %  [15,  9,  8,  7, 11]
    %  [10,  5, -1,  4,  6]
    %  [-1,  3,  2,  1, -1]
    %
    sensor_layout = [...
         -1, 14, 13, 12, -1;
         15,  9,  8,  7, 11;
         10,  5, -1,  4,  6;
         -1,  3,  2,  1, -1];

    % Figure 1: TF plots for Y-channels
    figure('Name', 'Time-Frequency Plots (Y-Channels)', 'Color', 'w', 'Position', [100, 100, 1200, 800]);
    plotTimeFrequencyWithLayout(avg_trials_y_bc, sf, TW, TS, meg_channel_names(y_idx), 'Time-Frequency Analysis (Y-Channels)', sensor_layout);

    % Figure 2: TF plots for Z-channels
    figure('Name', 'Time-Frequency Plots (Z-Channels)', 'Color', 'w', 'Position', [100, 100, 1200, 800]);
    plotTimeFrequencyWithLayout(avg_trials_z_bc, sf, TW, TS, meg_channel_names(z_idx), 'Time-Frequency Analysis (Z-Channels)', sensor_layout);


function plotTimeFrequencyWithLayout(data, fs, TW, TS, channel_names, title_str, sensor_layout)
    % data: [num_channels x num_timepoints]
    % sensor_layout: matrix (rows x cols) with sensor numbers (or -1 for an empty slot)
    
    num_channels = size(data, 1);
    num_timepoints = size(data, 2);
    nfft = TW;
    winfun = blackmanharris(TW);
    TWs = TW:TS:num_timepoints;
    freq = (0:(nfft/2)-1) * (fs / nfft);
    freq_limit = 70;          % Limit frequency to 100 Hz
    freq_idx = freq <= freq_limit;
    freq = freq(freq_idx);
    num_TW = numel(TWs);
    
    % Preallocate power spectrum array: [time windows x frequencies x channels]
    ps_on = zeros(num_TW, numel(freq), num_channels);
    
    % Compute time-frequency representation for each channel
    for ch = 1:num_channels
        disp(['Processing Channel ' num2str(ch) ' / ' num2str(num_channels)]);
        for tstep = 1:num_TW
            segment_start = TWs(tstep) - TW + 1;
            segment_end = TWs(tstep);
            if segment_start > 0 && segment_end <= num_timepoints
                dat = data(ch, segment_start:segment_end);
                dat = dat .* winfun';
                f_val = abs(fft(dat, nfft));
                ps_on(tstep, :, ch) = f_val(freq_idx);
            end
        end
    end
    
    % Baseline correction (using first 3 time windows)
    baseline_ps = mean(ps_on(1:3, :, :), 1);
    baseline_corrected_ps = ps_on ./ baseline_ps;
    log_baseline_corrected_ps = log10(max(baseline_corrected_ps, eps)); % Avoid log(0)
    time_freq_axis = (TWs - TW/2) / fs * 1000; % Time vector in ms
    
    % Determine subplot grid dimensions from sensor_layout
    [num_rows, num_cols] = size(sensor_layout);
    
    % Loop over the layout positions and plot if a sensor is assigned
    for row = 1:num_rows
        for col = 1:num_cols
            sensor_number = sensor_layout(row, col);
            if sensor_number == -1 || sensor_number > num_channels
                continue; % Skip empty or invalid positions
            end
            % Determine subplot index (row-major order)
            subplot_index = (row - 1) * num_cols + col;
            subplot(num_rows, num_cols, subplot_index);
            imagesc(time_freq_axis, freq, squeeze(log_baseline_corrected_ps(:, :, sensor_number))');
            axis tight;
            set(gca, 'YDir', 'normal');
            xlabel('Time (ms)');
            ylabel('Frequency (Hz)');
            title(channel_names{sensor_number});
            colorbar;
            clim([-2 2]); 
            hold on;
            % Mark trigger location (501 ms)
            plot([501 501], [min(freq) max(freq)], 'r', 'LineWidth', 1);
            hold off;
        end
    end
    
    sgtitle(title_str);
    % Add orientation annotations: Nose (top center), Left (left side), Right (right side)
    annotation('textbox', [0.45, 0.90, 0.1, 0.05], 'String', 'Nose', 'EdgeColor', 'none', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    annotation('textbox', [0.05, 0.45, 0.1, 0.05], 'String', 'Left', 'EdgeColor', 'none', 'FontSize', 14, 'FontWeight', 'bold', 'Rotation', 90, 'HorizontalAlignment', 'center');
    annotation('textbox', [0.90, 0.55, 0.1, 0.05], 'String', 'Right', 'EdgeColor', 'none', 'FontSize', 14, 'FontWeight', 'bold', 'Rotation', -90, 'HorizontalAlignment', 'center');
end
