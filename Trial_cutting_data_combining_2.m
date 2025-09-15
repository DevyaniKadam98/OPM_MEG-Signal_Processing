function Trial_cutting_data_combining(run_files, output_file)
    % Function to process OPM-MEG data for multiple runs, cut trials around triggers,
    % and combine them into a single dataset.

    % Parameters for trial cutting
    pre_trigger_time = 0.5; % Seconds before the trigger
    post_trigger_time = 1; % Seconds after the trigger

    % Initialize variables for concatenation
    combined_trials = [];
    combined_trial_time_axis = [];
    Fs = [];
    meg_channel_names = [];

    % Loop through each file
    for run_idx = 1:numel(run_files)
        % Load the raw data for the current run
        disp(['Processing file: ', run_files{run_idx}]);
        data = load(run_files{run_idx});

        % Extract CAR-corrected MEG signals and metadata
        meg_signals_car = data.meg_signals_car; % CAR-corrected MEG data
        trigger_pulse = data.trigger_pulse; % Third EMG channel is 'trigger_pulse'
        Fs_run = data.Fs; % Sampling frequency

        % Ensure consistency of metadata across runs
        if isempty(Fs)
            Fs = Fs_run;
        elseif Fs ~= Fs_run
            error('Sampling frequency mismatch in file: %s', run_files{run_idx});
        end

        if isempty(meg_channel_names)
            meg_channel_names = data.meg_channel_names;
        elseif ~isequal(meg_channel_names, data.meg_channel_names)
            error('Channel names mismatch in file: %s', run_files{run_idx});
        end

        % Convert times to samples
        pre_samples = round(pre_trigger_time * Fs);
        post_samples = round(post_trigger_time * Fs);
        trial_length = pre_samples + post_samples + 1; % Total length of each trial in samples

        % Find indices where the trigger pulse is 1
        trigger_indices = find(trigger_pulse == 1);

        % Initialize array to store trials for the current run
        num_trials = length(trigger_indices);
        num_channels = size(meg_signals_car, 1);
        trials_run = zeros(num_channels, trial_length, num_trials);

        % Cut trials around triggers
        for i = 1:num_trials
            start_idx = max(trigger_indices(i) - pre_samples, 1); % Ensure within bounds
            end_idx = min(trigger_indices(i) + post_samples, size(meg_signals_car, 2)); % Ensure within bounds

            % Check if trial length matches expected length
            if end_idx - start_idx + 1 < trial_length
                warning('Skipping trial %d in run %d: insufficient data at edges.', i, run_idx);
                continue;
            end

            % Extract trial
            trials_run(:, :, i) = meg_signals_car(:, start_idx:end_idx);
        end

        % Concatenate trials across runs
        combined_trials = cat(3, combined_trials, trials_run);

        % Define trial time axis (consistent across runs)
        if isempty(combined_trial_time_axis)
            combined_trial_time_axis = linspace(-pre_trigger_time, post_trigger_time, trial_length) * 1000; % Time in ms
        end
    end

    % Save the combined data
    save(output_file, 'combined_trials', 'meg_channel_names', 'combined_trial_time_axis', 'Fs', '-v7.3');
    disp(['Combined trials saved to "', output_file, '".']);
end

%  Usage:
run_files = {
    'CAR_Corrected_MEG_with_Processed_EMG_093_1.mat', ...
    'CAR_Corrected_MEG_with_Processed_EMG_093_2.mat', ...
    'CAR_Corrected_MEG_with_Processed_EMG_093_3.mat'
};
output_file = 'Cut_Trials_MEG_combined_093.mat';

% Process and combine trials
process_and_combine_trials(run_files, output_file);
