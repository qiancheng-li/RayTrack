function [nearest_positions, nearest_indices] = find_nearest_positions( ...
    query_positions, candidate_x, candidate_y)
%FIND_NEAREST_POSITIONS Match each query to its nearest candidate position.

    num_queries = size(query_positions, 1);
    nearest_positions = zeros(num_queries, 2);
    nearest_indices = zeros(num_queries, 1);

    for query_id = 1:num_queries
        x_difference = candidate_x - query_positions(query_id, 1);
        y_difference = candidate_y - query_positions(query_id, 2);
        candidate_distances = hypot(x_difference, y_difference);

        [~, nearest_index] = min(candidate_distances);
        nearest_positions(query_id, :) = [ ...
            candidate_x(nearest_index), candidate_y(nearest_index)];
        nearest_indices(query_id) = nearest_index;
    end
end
