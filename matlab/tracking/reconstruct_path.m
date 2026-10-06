function path = reconstruct_path(candidate_grid_indices, start_grid_index, ...
    end_grid_index, measured_rssi, recursion_depth, power_map, ...
    longitude_step, latitude_step, min_longitude, min_latitude)
%RECONSTRUCT_PATH Estimate a path by recursively selecting RSSI pivots.
%
% Each candidate pivot defines a two-line path between the known endpoints.
% The candidate whose simulated RSSI best matches the measurements is used
% to split the path, and the two halves are refined recursively.

    num_path_samples = numel(measured_rssi);

    if num_path_samples == 1
        path = start_grid_index;
        return;
    elseif num_path_samples == 2
        path = [start_grid_index; end_grid_index];
        return;
    end

    candidate_loss_map = inf(size(power_map));
    start_latlon = grid_index_to_latlon(start_grid_index, longitude_step, ...
        latitude_step, min_longitude, min_latitude);
    end_latlon = grid_index_to_latlon(end_grid_index, longitude_step, ...
        latitude_step, min_longitude, min_latitude);

    for candidate_id = 1:size(candidate_grid_indices, 1)
        candidate_grid_index = candidate_grid_indices(candidate_id, :);
        candidate_latlon = grid_index_to_latlon(candidate_grid_index, ...
            longitude_step, latitude_step, min_longitude, min_latitude);

        distance_from_start = geographic_distance_meters( ...
            start_latlon(1), start_latlon(2), ...
            candidate_latlon(1), candidate_latlon(2));
        distance_to_end = geographic_distance_meters( ...
            end_latlon(1), end_latlon(2), ...
            candidate_latlon(1), candidate_latlon(2));

        [first_half_length, second_half_length] = split_lengths( ...
            num_path_samples, distance_from_start, distance_to_end);

        candidate_path = interpolate_through_pivot(start_grid_index, ...
            candidate_grid_index, end_grid_index, ...
            first_half_length, second_half_length);
        candidate_path = round(candidate_path);

        simulated_rssi = sample_power_map(candidate_path, power_map);
        speed_penalty = calculate_speed_penalty( ...
            distance_from_start + distance_to_end, num_path_samples);
        candidate_loss = sum(abs(simulated_rssi - measured_rssi)) * ...
            speed_penalty;

        candidate_row = candidate_grid_index(1);
        candidate_column = candidate_grid_index(2);
        candidate_loss_map(candidate_row, candidate_column) = candidate_loss;
    end

    score_map = 1 ./ candidate_loss_map;
    score_map(~isfinite(score_map)) = 0;
    score_map = gaussian_smooth(score_map, 6);

    [~, best_linear_index] = max(score_map(:));
    [pivot_row, pivot_column] = ind2sub(size(score_map), best_linear_index);
    pivot_grid_index = [pivot_row, pivot_column];
    pivot_latlon = grid_index_to_latlon(pivot_grid_index, longitude_step, ...
        latitude_step, min_longitude, min_latitude);

    distance_from_start = geographic_distance_meters( ...
        start_latlon(1), start_latlon(2), ...
        pivot_latlon(1), pivot_latlon(2));
    distance_to_end = geographic_distance_meters( ...
        end_latlon(1), end_latlon(2), ...
        pivot_latlon(1), pivot_latlon(2));
    [first_half_length, second_half_length] = split_lengths( ...
        num_path_samples, distance_from_start, distance_to_end);

    if recursion_depth == 1
        path = interpolate_through_pivot(start_grid_index, pivot_grid_index, ...
            end_grid_index, first_half_length, second_half_length);
        return;
    end

    first_path = reconstruct_path(candidate_grid_indices, start_grid_index, ...
        pivot_grid_index, measured_rssi(1:first_half_length), ...
        recursion_depth - 1, power_map, longitude_step, latitude_step, ...
        min_longitude, min_latitude);
    second_path = reconstruct_path(candidate_grid_indices, pivot_grid_index, ...
        end_grid_index, measured_rssi(first_half_length:end), ...
        recursion_depth - 1, power_map, longitude_step, latitude_step, ...
        min_longitude, min_latitude);

    path = [first_path(1:end-1, :); second_path];
end

function [first_half_length, second_half_length] = split_lengths( ...
    total_length, distance_from_start, distance_to_end)
    total_distance = distance_from_start + distance_to_end;
    if total_distance == 0
        first_half_length = floor(total_length / 2);
    else
        first_half_length = floor( ...
            total_length * distance_from_start / total_distance);
    end
    first_half_length = max(first_half_length, 1);
    second_half_length = total_length - first_half_length + 1;
end

function path = interpolate_through_pivot( ...
    start_index, pivot_index, end_index, ...
    first_half_length, second_half_length)
    first_rows = linspace(start_index(1), pivot_index(1), first_half_length);
    second_rows = linspace(pivot_index(1), end_index(1), second_half_length);
    first_columns = linspace(start_index(2), pivot_index(2), first_half_length);
    second_columns = linspace(pivot_index(2), end_index(2), second_half_length);

    path = [ ...
        [first_rows, second_rows(2:end)].', ...
        [first_columns, second_columns(2:end)].' ...
    ];
end

function penalty = calculate_speed_penalty(path_distance, num_path_samples)
    standard_speed = 0.2624;
    estimated_speed = path_distance / max(num_path_samples - 1, 1);

    if estimated_speed == 0
        penalty = inf;
        return;
    end

    relative_speed = estimated_speed / standard_speed;
    penalty = max(relative_speed, 1 / relative_speed);
end

function smoothed_values = gaussian_smooth(values, sigma)
    kernel_radius = ceil(3 * sigma);
    [kernel_x, kernel_y] = meshgrid( ...
        -kernel_radius:kernel_radius, -kernel_radius:kernel_radius);
    gaussian_kernel = exp( ...
        -(kernel_x .^ 2 + kernel_y .^ 2) / (2 * sigma ^ 2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel, "all");
    smoothed_values = conv2(values, gaussian_kernel, "same");
end
