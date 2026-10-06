function [peak_values, peak_rows, peak_columns] = find_2d_peaks(data, varargin)
%FIND_2D_PEAKS Find strict local maxima in a two-dimensional array.
%
% Supported name-value options:
%   MinPeakHeight   Minimum accepted peak value. Default: -Inf.
%   Threshold       Minimum difference from every neighbor. Default: 0.
%   MinPeakDistance Minimum Euclidean separation between retained peaks.

    validateattributes(data, {"numeric"}, {"2d", "nonempty"});

    parser = inputParser;
    addParameter(parser, "MinPeakHeight", -inf, ...
        @(value) isnumeric(value) && isscalar(value));
    addParameter(parser, "Threshold", 0, ...
        @(value) isnumeric(value) && isscalar(value) && value >= 0);
    addParameter(parser, "MinPeakDistance", 0, ...
        @(value) isnumeric(value) && isscalar(value) && value >= 0);
    parse(parser, varargin{:});
    options = parser.Results;

    finite_data = double(data);
    finite_data(~isfinite(finite_data)) = -inf;
    neighbor_maximum = -inf(size(finite_data));

    [num_rows, num_columns] = size(finite_data);
    for row_offset = -1:1
        for column_offset = -1:1
            if row_offset == 0 && column_offset == 0
                continue;
            end

            source_rows = max(1, 1-row_offset): ...
                min(num_rows, num_rows-row_offset);
            source_columns = max(1, 1-column_offset): ...
                min(num_columns, num_columns-column_offset);
            target_rows = source_rows + row_offset;
            target_columns = source_columns + column_offset;

            shifted_neighbor = -inf(size(finite_data));
            shifted_neighbor(target_rows, target_columns) = ...
                finite_data(source_rows, source_columns);
            neighbor_maximum = max(neighbor_maximum, shifted_neighbor);
        end
    end

    is_peak = isfinite(finite_data) & ...
        finite_data > neighbor_maximum & ...
        finite_data >= options.MinPeakHeight & ...
        finite_data - neighbor_maximum >= options.Threshold;

    [peak_rows, peak_columns] = find(is_peak);
    peak_values = data(is_peak);

    [peak_values, sort_order] = sort(peak_values, "descend");
    peak_rows = peak_rows(sort_order);
    peak_columns = peak_columns(sort_order);

    if options.MinPeakDistance > 0
        keep_peak = false(size(peak_values));
        selected_rows = zeros(0, 1);
        selected_columns = zeros(0, 1);

        for peak_id = 1:numel(peak_values)
            distances = hypot( ...
                peak_rows(peak_id) - selected_rows, ...
                peak_columns(peak_id) - selected_columns);
            if isempty(distances) || all(distances >= options.MinPeakDistance)
                keep_peak(peak_id) = true;
                selected_rows(end+1, 1) = peak_rows(peak_id); %#ok<AGROW>
                selected_columns(end+1, 1) = peak_columns(peak_id); %#ok<AGROW>
            end
        end

        peak_values = peak_values(keep_peak);
        peak_rows = peak_rows(keep_peak);
        peak_columns = peak_columns(keep_peak);
    end
end
