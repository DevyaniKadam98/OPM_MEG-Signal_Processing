function plot_avg_trials_with_trigger()
    % Load the data file
    data = load('Cut_Trials_MEG_combined_093.mat');
    
    % Extract relevant data
    trials = data.combined_trials; % Trials data: [sensors x time points x trials]
    trial_time_axis = data.combined_trial_time_axis; % Time axis for trials in ms
    meg_channel_names = data.meg_channel_names; % Sensor names

    % Identify Y-axis and Z-axis sensors
    y_sensors_idx = find(contains(meg_channel_names, '_y')); % Indices of Y-axis sensors
    z_sensors_idx = find(contains(meg_channel_names, '_z')); % Indices of Z-axis sensors

    % Average the trials for each sensor
    avg_trials = mean(trials, 3) * 1e3; % Average across the 3rd dimension (convert nT to pT)

    % Define trigger position in time axis
    trigger_position = 0; % Trigger is at time 0 (relative position in the trial)
    trigger_index = find(trial_time_axis >= trigger_position, 1, 'first'); % Index of the trigger

    % Determine global y-axis limits for consistent scaling
    y_data = avg_trials(y_sensors_idx, :); % Y-axis sensor data
    z_data = avg_trials(z_sensors_idx, :); % Z-axis sensor data
    global_min = min([y_data(:); z_data(:)]);
    global_max = max([y_data(:); z_data(:)]);

    % Plot all Y-axis sensors in one subplot
    figure('Name', 'Averaged Trials with Trigger Position', 'Color', 'w', 'Position', [100, 100, 1200, 600]);

    subplot(2, 1, 1);
    hold on;
    for i = 1:numel(y_sensors_idx)
        sensor_idx = y_sensors_idx(i);
        plot(trial_time_axis, avg_trials(sensor_idx, :), 'DisplayName', meg_channel_names{sensor_idx}, 'LineWidth', 1.5);
    end
    xline(trial_time_axis(trigger_index), 'r--', 'LineWidth', 1.5); % Trigger position
    title('Averaged Trials for Y-Axis Sensors');
    xlabel('Time (ms)');
    ylabel('Field (pT)'); % Field in picotesla
    ylim([-20, 20]); % Apply consistent y-scale
    yticks(-20:5:20);
    legend('show');
    grid on;
    hold off;

    % Plot all Z-axis sensors in another subplot
    subplot(2, 1, 2);
    hold on;
    for i = 1:numel(z_sensors_idx)
        sensor_idx = z_sensors_idx(i);
        plot(trial_time_axis, avg_trials(sensor_idx, :), 'DisplayName', meg_channel_names{sensor_idx}, 'LineWidth', 1.5);
    end
    xline(trial_time_axis(trigger_index), 'r--', 'LineWidth', 1.5); % Trigger position
    title('Averaged Trials for Z-Axis Sensors');
    xlabel('Time (ms)');
    ylabel('Field (pT)'); % Field in picotesla
    ylim([-30, 30]); % Apply consistent y-scale
    yticks(-30:5:30);
    legend('show');
    grid on;
    hold off;

    sgtitle('Averaged Trials with Trigger Position (Subject093 All runs)');
end
