function power_values = sample_power_map(grid_indices, power_map)
%SAMPLE_POWER_MAP Return radio-map values at [row, column] grid indices.

    row_indices = grid_indices(:, 1);
    column_indices = grid_indices(:, 2);

    if any(row_indices < 1 | row_indices > size(power_map, 1) | ...
            column_indices < 1 | column_indices > size(power_map, 2))
        error("sample_power_map:IndexOutOfBounds", ...
            "At least one grid index lies outside the power map.");
    end

    linear_indices = sub2ind(size(power_map), row_indices, column_indices);
    power_values = power_map(linear_indices);
end
