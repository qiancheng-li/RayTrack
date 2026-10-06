function [median_error, error_std, position_errors, segment_times] = ...
    evaluate_tracking_path(data_file, samples_per_segment, recursion_depth, base_station_name)
%EVALUATE_TRACKING_PATH Reconstruct and evaluate a trajectory from RSSI.
%
% The measured trajectory is divided into overlapping segments. Each
% segment has a known start and end position, and RECONSTRUCT_PATH estimates
% the intermediate points by matching measured and simulated RSSI.

    arguments
        data_file
        samples_per_segment (1, 1) double {mustBeInteger, mustBePositive}
        recursion_depth (1, 1) double {mustBeInteger, mustBePositive}
        base_station_name
    end

    coverage_map_dir = fullfile(".", "Coverage Map");
    power_grid_file = fullfile(coverage_map_dir, ...
        string(base_station_name) + "_power_grid.mat");
    power_map_file = fullfile(coverage_map_dir, ...
        "power_square_" + string(base_station_name) + ".mat");

    % Linear calibration from measured RSSI to simulated RSSI.
    calibration_slope = 0.56745;
    calibration_offset = 17.8535;

    measurements = readmatrix(data_file);
    measured_rssi = measurements(:, 3) * calibration_slope + calibration_offset;

    power_grid = importdata(power_grid_file);
    power_map = importdata(power_map_file);

    map_latitudes = power_grid(:, 1);
    map_longitudes = power_grid(:, 2);
    min_latitude = min(map_latitudes);
    min_longitude = min(map_longitudes);

    unique_latitudes = unique(map_latitudes, "stable");
    unique_longitudes = unique(map_longitudes, "stable");
    latitude_step = abs(unique_latitudes(1) - unique_latitudes(2));
    longitude_step = abs(unique_longitudes(1) - unique_longitudes(2));

    map_positions = power_grid(:, 1:2);
    map_grid_indices = round(latlon_to_grid_index(map_positions, longitude_step, ...
        latitude_step, min_longitude, min_latitude));

    num_samples = numel(measured_rssi);
    num_full_segments = floor(num_samples / samples_per_segment);
    estimated_path = zeros(num_samples, 2);
    segment_times = zeros(num_full_segments, 1);

    for segment_id = 1:num_full_segments
        fprintf("Tracking segment %d of %d\n", segment_id, num_full_segments);

        segment_start = (segment_id - 1) * samples_per_segment + 1;
        next_segment_start = segment_id * samples_per_segment + 1;

        start_latlon = measurements(segment_start, 1:2);
        start_grid_index = latlon_to_grid_index(start_latlon, longitude_step, ...
            latitude_step, min_longitude, min_latitude);

        if next_segment_start <= num_samples
            segment_end = next_segment_start;
        else
            segment_end = segment_id * samples_per_segment;
        end

        end_latlon = measurements(segment_end, 1:2);
        end_grid_index = latlon_to_grid_index(end_latlon, longitude_step, ...
            latitude_step, min_longitude, min_latitude);

        candidate_rows = find_candidate_points(start_latlon, end_latlon, ...
            map_positions);
        candidate_grid_indices = round(latlon_to_grid_index( ...
            map_positions(candidate_rows, :), longitude_step, latitude_step, ...
            min_longitude, min_latitude));

        tic;
        segment_path = reconstruct_path(candidate_grid_indices, ...
            start_grid_index, end_grid_index, ...
            measured_rssi(segment_start:segment_end), recursion_depth, ...
            power_map, longitude_step, latitude_step, ...
            min_longitude, min_latitude);
        segment_times(segment_id) = toc;

        if next_segment_start <= num_samples
            estimated_path(segment_start:segment_end-1, :) = ...
                segment_path(1:end-1, :);
        else
            estimated_path(segment_start:segment_end, :) = segment_path;
        end
    end

    % Process samples left after the final complete segment.
    remaining_start = num_full_segments * samples_per_segment + 1;
    if remaining_start <= num_samples
        start_latlon = measurements(remaining_start, 1:2);
        end_latlon = measurements(end, 1:2);

        start_grid_index = latlon_to_grid_index(start_latlon, longitude_step, ...
            latitude_step, min_longitude, min_latitude);
        end_grid_index = latlon_to_grid_index(end_latlon, longitude_step, ...
            latitude_step, min_longitude, min_latitude);

        candidate_rows = find_candidate_points(start_latlon, end_latlon, ...
            map_positions);
        candidate_grid_indices = map_grid_indices(candidate_rows, :);

        estimated_path(remaining_start:end, :) = reconstruct_path( ...
            candidate_grid_indices, start_grid_index, end_grid_index, ...
            measured_rssi(remaining_start:end), recursion_depth, power_map, ...
            longitude_step, latitude_step, min_longitude, min_latitude);
    end

    estimated_latlon = grid_index_to_latlon(estimated_path, longitude_step, ...
        latitude_step, min_longitude, min_latitude);

    position_errors = geographic_distance_meters(estimated_latlon(:, 1), ...
        estimated_latlon(:, 2), measurements(:, 1), measurements(:, 2));

    median_error = median(position_errors);
    error_std = std(position_errors);
end
