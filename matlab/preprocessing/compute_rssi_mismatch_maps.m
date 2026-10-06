function mismatch_maps = compute_rssi_mismatch_maps( ...
    simulated_rssi_map, measured_rssi)
%COMPUTE_RSSI_MISMATCH_MAPS Compare measured RSSI with every radio-map cell.
%
% Inputs
%   simulated_rssi_map : [rows, columns, frequencies]
%   measured_rssi      : [locations, frequencies]
%
% Output
%   mismatch_maps      : [rows, columns, locations], containing the sum of
%                        absolute RSSI differences across frequencies.

    validateattributes(simulated_rssi_map, {"numeric"}, {"nonempty"});
    if ndims(simulated_rssi_map) ~= 3
        error("compute_rssi_mismatch_maps:InvalidRadioMapShape", ...
            "simulated_rssi_map must have dimensions [rows, columns, frequencies].");
    end
    validateattributes(measured_rssi, {"numeric"}, ...
        {"nonempty", "2d", "finite"});

    num_frequencies = size(simulated_rssi_map, 3);
    if size(measured_rssi, 2) ~= num_frequencies
        error("compute_rssi_mismatch_maps:FrequencyCountMismatch", ...
            ["The radio map has %d frequency slices, but the measured " ...
             "data has %d frequencies."], ...
            num_frequencies, size(measured_rssi, 2));
    end

    map_rows = size(simulated_rssi_map, 1);
    map_columns = size(simulated_rssi_map, 2);
    num_locations = size(measured_rssi, 1);

    simulated_rssi_map = single(simulated_rssi_map);
    measured_rssi = single(measured_rssi);
    mismatch_maps = inf( ...
        map_rows, map_columns, num_locations, "single");

    for location_id = 1:num_locations
        measured_vector = reshape( ...
            measured_rssi(location_id, :), 1, 1, num_frequencies);
        mismatch_maps(:, :, location_id) = sum( ...
            abs(simulated_rssi_map - measured_vector), 3);
    end
end
