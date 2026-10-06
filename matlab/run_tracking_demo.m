%% Minimal tracking demo using only data included in this repository
clear;
clc;
close all;

%% Resolve repository paths and choose a short demonstration segment
matlab_dir = fileparts(mfilename("fullpath"));
repo_root = fileparts(matlab_dir);
addpath(genpath(matlab_dir));

radio_map_file = fullfile(repo_root, "data", "algorithm_inputs", ...
    "lab1_radio_map_points.mat");
measurement_file = fullfile(repo_root, "data", ...
    "experimental_standardized", "Station-Lab1", "lab1 data", ...
    "11-14-1-interpolated.csv");

demo_sample_count = 21;
recursion_depth = 2;
% Match the calibration used by the full tracking evaluator.
calibration_slope = 0.56745;
calibration_offset = 17.8535;

%% Load the included radio-map points and build a regular grid
radio_map_data = load(radio_map_file);
radio_map_points = radio_map_data.power_map_grid;

map_positions = radio_map_points(:, 1:2);
map_rssi = radio_map_points(:, 3);

unique_latitudes = sort(unique(map_positions(:, 1)));
unique_longitudes = sort(unique(map_positions(:, 2)));
latitude_step = median(diff(unique_latitudes));
longitude_step = median(diff(unique_longitudes));
min_latitude = min(unique_latitudes);
min_longitude = min(unique_longitudes);

valid_grid_indices = round(latlon_to_grid_index(map_positions, longitude_step, ...
    latitude_step, min_longitude, min_latitude));
num_rows = max(valid_grid_indices(:, 1));
num_columns = max(valid_grid_indices(:, 2));

% Fill cells outside the sparse valid set with the nearest available RSSI,
% allowing every interpolated candidate path to be sampled.
power_map = nan(num_rows, num_columns);
valid_linear_indices = sub2ind(size(power_map), ...
    valid_grid_indices(:, 1), valid_grid_indices(:, 2));
power_map(valid_linear_indices) = map_rssi;
power_map = fillmissing(power_map, "nearest", 1, "EndValues", "nearest");
power_map = fillmissing(power_map, "nearest", 2, "EndValues", "nearest");

%% Load a short measured trajectory
measurements = readmatrix(measurement_file);
measurements = measurements( ...
    1:min(demo_sample_count, size(measurements, 1)), :);
measured_rssi = measurements(:, 3) * calibration_slope + calibration_offset;

start_latlon = measurements(1, 1:2);
end_latlon = measurements(end, 1:2);
start_grid_index = latlon_to_grid_index(start_latlon, longitude_step, latitude_step, ...
    min_longitude, min_latitude);
end_grid_index = latlon_to_grid_index(end_latlon, longitude_step, latitude_step, ...
    min_longitude, min_latitude);

candidate_rows = find_candidate_points(start_latlon, end_latlon, ...
    map_positions);
if isempty(candidate_rows)
    error("No candidate radio-map points were found for the demo segment.");
end
candidate_grid_indices = valid_grid_indices(candidate_rows, :);

%% Reconstruct and evaluate the trajectory
estimated_grid_path = reconstruct_path(candidate_grid_indices, ...
    start_grid_index, end_grid_index, measured_rssi, recursion_depth, ...
    power_map, longitude_step, latitude_step, ...
    min_longitude, min_latitude);

estimated_latlon = grid_index_to_latlon(estimated_grid_path, longitude_step, ...
    latitude_step, min_longitude, min_latitude);
position_errors = geographic_distance_meters(estimated_latlon(:, 1), ...
    estimated_latlon(:, 2), measurements(:, 1), measurements(:, 2));

fprintf("Tracking demo complete.\n");
fprintf("Samples: %d, candidate points: %d\n", ...
    size(measurements, 1), numel(candidate_rows));
fprintf("Median error: %.2f m\n", median(position_errors));
fprintf("Mean error: %.2f m\n", mean(position_errors));

figure("Color", "w");
plot(measurements(:, 2), measurements(:, 1), ...
    "o-", "LineWidth", 1.5, "DisplayName", "Measured trajectory");
hold on;
plot(estimated_latlon(:, 2), estimated_latlon(:, 1), ...
    "x-", "LineWidth", 1.5, "DisplayName", "Estimated trajectory");
axis equal;
grid on;
xlabel("Longitude");
ylabel("Latitude");
title("Minimal RSSI tracking demo");
legend("Location", "best");
