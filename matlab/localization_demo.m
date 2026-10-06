%% Single-base-station localization demo for Lab 1
clear;
clc;
close all;

matlab_dir = fileparts(mfilename("fullpath"));
repo_root = fileparts(matlab_dir);
addpath(fullfile(matlab_dir, "helpers"));

%% Configuration
lab1_file_id = 1;       % freq_map1.mat and DISTMASTER27freq1.mat
num_evaluations = 100;
frequency_slice = 14;   % 915 MHz slice in the original radio map

rssi_offset = 0;
motion_offset = 0;
path_score_offset = 0;
use_path_constraint = true;
recent_path_length = 10;

origin_latitude = 47.655185;
origin_longitude = -122.307035;
grid_row_offset = 692.51605;
grid_column_offset = 409.3338;

window_latitude = 47.653750;
window_longitude = -122.306494;
vertical_distance = geographic_distance_meters( ...
    window_latitude, origin_longitude, ...
    origin_latitude, origin_longitude);
horizontal_distance = geographic_distance_meters( ...
    origin_latitude, window_longitude, ...
    origin_latitude, origin_longitude);
latitude_step = ...
    (window_latitude - origin_latitude) / vertical_distance;
longitude_step = ...
    (window_longitude - origin_longitude) / horizontal_distance;

%% Load inputs
latitude_data = load(fullfile(repo_root, "data", "algorithm_inputs", ...
    "fixed_location_latitudes.mat"));
longitude_data = load(fullfile(repo_root, "data", "algorithm_inputs", ...
    "fixed_location_longitudes.mat"));
mismatch_data = load(fullfile(repo_root, "freqDistance", ...
    sprintf("DISTMASTER27freq%d.mat", lab1_file_id)));
last_location_data = load(fullfile( ...
    repo_root, "freqDistance", "DISTMASTERlastloc.mat"));
recent_path_data = load(fullfile(repo_root, "data", "algorithm_inputs", ...
    "lab1_last_path_rssi_915mhz.mat"));
radio_map_data = load(fullfile( ...
    repo_root, sprintf("freq_map%d.mat", lab1_file_id)));

target_latitudes = latitude_data.datalats;
target_longitudes = longitude_data.datalons;
rssi_mismatch_maps = mismatch_data.DISTMASTER;
last_location_maps = last_location_data.lastloc_dis;
recent_path_rssi = recent_path_data.last_paths;
lab1_power_map = radio_map_data.freqMap(:, :, frequency_slice);
map_size = size(lab1_power_map);

num_evaluations = min(num_evaluations, numel(target_latitudes));
if size(rssi_mismatch_maps, 3) < num_evaluations
    error("localization_demo:InsufficientMismatchMaps", ...
        "RSSI mismatch data contains fewer locations than requested.");
end
if size(last_location_maps, 3) < num_evaluations
    error("localization_demo:InsufficientMotionMaps", ...
        "Last-location data contains fewer locations than requested.");
end
if size(recent_path_rssi, 1) < num_evaluations || ...
        size(recent_path_rssi, 2) < recent_path_length + 1
    error("localization_demo:InsufficientPathData", ...
        "Recent-path data does not match the requested evaluation size.");
end

%% Localize each fixed position
localization_errors = zeros(num_evaluations, 1);
predicted_grid_indices = zeros(num_evaluations, 2);
actual_grid_indices = zeros(num_evaluations, 2);

for location_id = 1:num_evaluations
    cost_map = ...
        (rssi_mismatch_maps(:, :, location_id) + rssi_offset) .* ...
        (last_location_maps(:, :, location_id) + motion_offset + 1e-12);
    score_map = smoothdata(1 ./ (cost_map + 1e-12), "sgolay", 25);

    peak_threshold = 0.1 * max(score_map(:));
    [peak_values, peak_rows, peak_columns] = find_2d_peaks( ...
        score_map, "MinPeakHeight", peak_threshold);

    if isempty(peak_values)
        [peak_values, linear_index] = max(score_map(:));
        [peak_rows, peak_columns] = ind2sub(map_size, linear_index);
    end

    if use_path_constraint
        path_scores = zeros(numel(peak_values), 1);
        measured_path = squeeze(recent_path_rssi(location_id, :, 1));

        for peak_id = 1:numel(peak_values)
            peak_row = peak_rows(peak_id);
            peak_column = peak_columns(peak_id);
            best_path_error = inf;

            for start_row = ...
                    peak_row-recent_path_length:peak_row+recent_path_length
                if start_row < 1 || start_row > map_size(1)
                    continue;
                end

                for start_column = ...
                        peak_column-recent_path_length: ...
                        peak_column+recent_path_length
                    if start_column < 1 || start_column > map_size(2)
                        continue;
                    end

                    start_latlon = sionna_grid_index_to_latlon( ...
                        [start_row, start_column], ...
                        longitude_step, latitude_step, ...
                        origin_longitude, origin_latitude);
                    peak_latlon = sionna_grid_index_to_latlon( ...
                        [peak_row, peak_column], ...
                        longitude_step, latitude_step, ...
                        origin_longitude, origin_latitude);
                    candidate_distance = geographic_distance_meters( ...
                        start_latlon(1), start_latlon(2), ...
                        peak_latlon(1), peak_latlon(2));

                    if abs(candidate_distance - recent_path_length) > 1
                        continue;
                    end

                    candidate_rows = round(linspace( ...
                        start_row, peak_row, recent_path_length + 1));
                    candidate_columns = round(linspace( ...
                        start_column, peak_column, recent_path_length + 1));
                    simulated_path = sample_power_map( ...
                        [candidate_rows.', candidate_columns.'], ...
                        lab1_power_map);
                    path_error = sum( ...
                        abs(simulated_path - measured_path), "omitnan");
                    best_path_error = min(best_path_error, path_error);
                end
            end

            path_scores(peak_id) = 1 / (best_path_error + 1e-12);
        end

        final_scores = ...
            peak_values .* (path_scores + path_score_offset);
    else
        final_scores = peak_values;
    end

    [~, best_peak_id] = max(final_scores);
    estimated_row = peak_rows(best_peak_id);
    estimated_column = peak_columns(best_peak_id);
    predicted_grid_indices(location_id, :) = ...
        [estimated_row, estimated_column];

    [actual_x, actual_y] = latlon_to_local_xy( ...
        target_latitudes(location_id), target_longitudes(location_id), ...
        origin_latitude, origin_longitude);
    actual_grid_index = [ ...
        actual_y + grid_row_offset, actual_x + grid_column_offset];
    actual_grid_indices(location_id, :) = actual_grid_index;
    localization_errors(location_id) = hypot( ...
        actual_grid_index(1) - estimated_row, ...
        actual_grid_index(2) - estimated_column);
end

%% Report and visualize results
fprintf("Single-BS Lab 1 localization result:\n");
fprintf("Median error: %.3f m\n", median(localization_errors));
fprintf("Mean error:   %.3f m\n", mean(localization_errors));
fprintf("Std error:    %.3f m\n", std(localization_errors));

figure("Color", "w");
plot_empirical_cdf(localization_errors);
grid on;
xlabel("Localization error (m)");
ylabel("CDF");
title("Single-BS Lab 1 localization error");

figure("Color", "w");
scatter(actual_grid_indices(:, 2), actual_grid_indices(:, 1), ...
    30, "red", "filled", "DisplayName", "Ground truth");
hold on;
scatter(predicted_grid_indices(:, 2), predicted_grid_indices(:, 1), ...
    30, "magenta", "filled", "DisplayName", "Prediction");
grid on;
axis equal;
legend("Location", "best");
xlabel("Grid column");
ylabel("Grid row");
title("Ground truth and predicted positions");
