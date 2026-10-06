function constraint_maps = compute_last_location_constraint_maps( ...
    map_size, previous_locations, expected_travel_distance, ...
    origin_latlon, index_offsets, valid_locations)
%COMPUTE_LAST_LOCATION_CONSTRAINT_MAPS Build motion-distance penalty maps.
%
% Each output cell is the absolute difference between its distance from the
% previous location and the expected travel distance.
%
% Inputs
%   map_size                 : [rows, columns]
%   previous_locations       : [samples, 2] latitude/longitude
%   expected_travel_distance : scalar distance in meters/grid cells
%   origin_latlon            : [origin_latitude, origin_longitude]
%   index_offsets            : [row_offset, column_offset]
%   valid_locations          : optional logical mask

    validateattributes(map_size, {"numeric"}, ...
        {"vector", "numel", 2, "integer", "positive"});
    validateattributes(previous_locations, {"numeric"}, ...
        {"2d", "ncols", 2});
    validateattributes(expected_travel_distance, {"numeric"}, ...
        {"scalar", "finite", "nonnegative"});

    num_locations = size(previous_locations, 1);
    if nargin < 6
        valid_locations = all(isfinite(previous_locations), 2);
    end
    valid_locations = logical(valid_locations(:));

    if numel(valid_locations) ~= num_locations
        error("compute_last_location_constraint_maps:MaskSizeMismatch", ...
            "valid_locations must contain one value per location.");
    end

    map_rows = map_size(1);
    map_columns = map_size(2);
    [row_grid, column_grid] = ndgrid(1:map_rows, 1:map_columns);
    constraint_maps = inf( ...
        map_rows, map_columns, num_locations, "single");

    for location_id = 1:num_locations
        if ~valid_locations(location_id)
            continue;
        end

        [x_meters, y_meters] = latlon_to_local_xy( ...
            previous_locations(location_id, 1), ...
            previous_locations(location_id, 2), ...
            origin_latlon(1), origin_latlon(2));

        previous_row = y_meters + index_offsets(1);
        previous_column = x_meters + index_offsets(2);
        distance_from_previous = hypot( ...
            row_grid - previous_row, column_grid - previous_column);

        constraint_maps(:, :, location_id) = single(abs( ...
            distance_from_previous - expected_travel_distance));
    end
end
