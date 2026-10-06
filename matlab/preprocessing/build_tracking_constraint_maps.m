%% Build last-location and recent-path constraints from Lab 1 tracking runs
clear;
clc;

%% Paths and matching settings
preprocessing_dir = fileparts(mfilename("fullpath"));
matlab_dir = fileparts(preprocessing_dir);
repo_root = fileparts(matlab_dir);
addpath(preprocessing_dir);
addpath(fullfile(matlab_dir, "helpers"));

algorithm_input_dir = fullfile(repo_root, "data", "algorithm_inputs");
tracking_dir = fullfile(repo_root, "data", ...
    "experimental_standardized", "Station-Lab1", "lab1 data");
output_dir = fullfile(repo_root, "freqDistance");

% The large Lab 1 multi-frequency map is intentionally distributed
% separately. Only its first two dimensions are needed here.
lab1_radio_map_file = fullfile(repo_root, "freq_map1.mat");

search_start_row = 61;
previous_offset_samples = 26;
recent_path_length = 10;
maximum_match_distance_meters = 2;

origin_latlon = [47.655185, -122.307035];
index_offsets = [692.51605, 409.3338]; % [row, column]

if ~isfolder(output_dir)
    mkdir(output_dir);
end
if ~isfile(lab1_radio_map_file)
    error("build_tracking_constraint_maps:MissingRadioMap", ...
        "Missing Lab 1 radio map: %s", lab1_radio_map_file);
end

latitude_data = load(fullfile( ...
    algorithm_input_dir, "fixed_location_latitudes.mat"));
longitude_data = load(fullfile( ...
    algorithm_input_dir, "fixed_location_longitudes.mat"));
fixed_locations = [latitude_data.datalats, longitude_data.datalons];
num_locations = size(fixed_locations, 1);

% Discover runs from the standardized folder instead of hard-coding dates.
track_files = dir(fullfile(tracking_dir, "*-interpolated.csv"));
if isempty(track_files)
    error("build_tracking_constraint_maps:NoTrackingFiles", ...
        "No interpolated tracking CSV files were found in %s.", tracking_dir);
end

track_names = string({track_files.name}).';
track_data = cell(numel(track_files), 1);
best_match_distance = inf(num_locations, 1);
best_track_id = zeros(num_locations, 1);
best_row_index = zeros(num_locations, 1);

%% Match every fixed position to the nearest tracking sample
for track_id = 1:numel(track_files)
    track_data{track_id} = readmatrix( ...
        fullfile(tracking_dir, track_files(track_id).name));
    current_track = track_data{track_id};

    if size(current_track, 1) < search_start_row
        warning("Skipping short tracking run: %s", track_files(track_id).name);
        continue;
    end

    search_rows = search_start_row:size(current_track, 1);
    [nearest_positions, nearest_local_indices] = find_nearest_positions( ...
        fixed_locations, current_track(search_rows, 1), ...
        current_track(search_rows, 2));
    nearest_global_indices = ...
        nearest_local_indices + search_start_row - 1;
    match_distances = geographic_distance_meters( ...
        nearest_positions(:, 1), nearest_positions(:, 2), ...
        fixed_locations(:, 1), fixed_locations(:, 2));

    improved = match_distances < best_match_distance;
    best_match_distance(improved) = match_distances(improved);
    best_track_id(improved) = track_id;
    best_row_index(improved) = nearest_global_indices(improved);
end

% The 2 m threshold is reported as quality metadata, not used to discard
% otherwise valid nearest-run matches.
valid_match = best_track_id > 0;
within_match_threshold = ...
    valid_match & best_match_distance <= maximum_match_distance_meters;
matched_locations = nan(num_locations, 2);
previous_locations = nan(num_locations, 2);
valid_previous_location = false(num_locations, 1);

for location_id = 1:num_locations
    if ~valid_match(location_id)
        continue;
    end

    current_track = track_data{best_track_id(location_id)};
    matched_row = best_row_index(location_id);
    previous_row = matched_row - previous_offset_samples;

    matched_locations(location_id, :) = current_track(matched_row, 1:2);
    if previous_row >= 1
        previous_locations(location_id, :) = ...
            current_track(previous_row, 1:2);
        valid_previous_location(location_id) = true;
    end
end

travel_distances = geographic_distance_meters( ...
    previous_locations(valid_previous_location, 1), ...
    previous_locations(valid_previous_location, 2), ...
    matched_locations(valid_previous_location, 1), ...
    matched_locations(valid_previous_location, 2));
if isempty(travel_distances)
    error("build_tracking_constraint_maps:NoValidPreviousLocations", ...
        "No matched positions have enough history for a previous location.");
end
expected_travel_distance = median(travel_distances, "omitnan");

%% Build the last-location constraint maps
radio_map_data = load(lab1_radio_map_file);
if ~isfield(radio_map_data, "freqMap")
    error("build_tracking_constraint_maps:MissingVariable", ...
        "%s does not contain a freqMap variable.", lab1_radio_map_file);
end
map_size = [size(radio_map_data.freqMap, 1), ...
    size(radio_map_data.freqMap, 2)];

lastloc_dis = compute_last_location_constraint_maps( ...
    map_size, previous_locations, expected_travel_distance, ...
    origin_latlon, index_offsets, valid_previous_location);
valid_lastloc = valid_previous_location;

last_location_file = fullfile(output_dir, "DISTMASTERlastloc.mat");
save(last_location_file, "lastloc_dis", "previous_locations", ...
    "matched_locations", "best_match_distance", "best_track_id", ...
    "best_row_index", "track_names", "valid_lastloc", ...
    "within_match_threshold", "maximum_match_distance_meters", ...
    "previous_offset_samples", "expected_travel_distance", "-v7.3");

%% Extract the recent Lab 1 RSSI path ending at each matched position
last_paths = nan( ...
    num_locations, recent_path_length + 1, 1, "single");
valid_lastpath = false(num_locations, 1);

for location_id = 1:num_locations
    if ~valid_match(location_id)
        continue;
    end

    current_track = track_data{best_track_id(location_id)};
    path_end_row = best_row_index(location_id);
    path_start_row = path_end_row - recent_path_length;

    if path_start_row < 1
        continue;
    end

    last_paths(location_id, :, 1) = single( ...
        current_track(path_start_row:path_end_row, 3).');
    valid_lastpath(location_id) = true;
end

last_path_length = recent_path_length;
last_path_file = fullfile( ...
    algorithm_input_dir, "lab1_last_path_rssi_915mhz.mat");
save(last_path_file, "last_paths", "valid_lastpath", ...
    "last_path_length", "best_track_id", "best_row_index", "track_names");

fprintf("Saved last-location constraints: %s\n", last_location_file);
fprintf("Saved recent-path RSSI: %s\n", last_path_file);
fprintf("Valid location matches: %d/%d\n", sum(valid_match), num_locations);
fprintf("Matches within %.1f m: %d/%d\n", ...
    maximum_match_distance_meters, sum(within_match_threshold), num_locations);
fprintf("Valid recent paths: %d/%d\n", sum(valid_lastpath), num_locations);
