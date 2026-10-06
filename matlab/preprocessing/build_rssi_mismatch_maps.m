%% Build multi-frequency RSSI mismatch maps for all four base stations
clear;
clc;

%% Paths and calibration
preprocessing_dir = fileparts(mfilename("fullpath"));
matlab_dir = fileparts(preprocessing_dir);
repo_root = fileparts(matlab_dir);
addpath(preprocessing_dir);

% The large freq_map0.mat ... freq_map3.mat files are intentionally not
% committed to Git. Place downloaded copies in the repository root or
% change radio_map_dir below.
radio_map_dir = repo_root;
output_dir = fullfile(repo_root, "freqDistance");
data_dir = fullfile(repo_root, "data", "experimental_standardized");

station_names = ["Roof", "Lab1", "Lab2", "Chem"];
calibration_slope = 1;
calibration_offset = 0;

if ~isfolder(output_dir)
    mkdir(output_dir);
end

%% Generate one mismatch volume per base station
for station_id = 0:numel(station_names)-1
    station_name = station_names(station_id + 1);
    radio_map_file = fullfile( ...
        radio_map_dir, sprintf("freq_map%d.mat", station_id));
    measurement_file = fullfile( ...
        data_dir, "Station-" + station_name, ...
        "localization_measurements.csv");

    if ~isfile(radio_map_file)
        error("build_rssi_mismatch_maps:MissingRadioMap", ...
            "Missing radio map: %s", radio_map_file);
    end
    if ~isfile(measurement_file)
        error("build_rssi_mismatch_maps:MissingMeasurements", ...
            "Missing localization measurements: %s", measurement_file);
    end

    radio_map_data = load(radio_map_file);
    if ~isfield(radio_map_data, "freqMap")
        error("build_rssi_mismatch_maps:MissingVariable", ...
            "%s does not contain a freqMap variable.", radio_map_file);
    end
    simulated_rssi_map = radio_map_data.freqMap;

    % Each CSV contains 100 positions with one row per frequency.
    measurements = readtable(measurement_file);
    measurements = sortrows(measurements, ["location_id", "sample_index"]);
    location_ids = unique(measurements.location_id, "stable");
    num_locations = numel(location_ids);
    num_frequencies = size(simulated_rssi_map, 3);

    measured_rssi = nan(num_locations, num_frequencies);
    frequency_hz = nan(1, num_frequencies);

    for location_index = 1:num_locations
        rows = measurements.location_id == location_ids(location_index);
        location_measurements = measurements(rows, :);

        if height(location_measurements) ~= num_frequencies
            error("build_rssi_mismatch_maps:IncompleteLocation", ...
                ["Station %s location %d has %d frequency samples; " ...
                 "the radio map has %d slices."], ...
                station_name, location_ids(location_index), ...
                height(location_measurements), num_frequencies);
        end

        measured_rssi(location_index, :) = ...
            location_measurements.rssi_dbm.';

        if location_index == 1
            frequency_hz = location_measurements.frequency_hz.';
        elseif any(location_measurements.frequency_hz.' ~= frequency_hz)
            error("build_rssi_mismatch_maps:FrequencyOrderMismatch", ...
                "Frequency ordering differs between localization positions.");
        end
    end

    measured_rssi = ...
        measured_rssi * calibration_slope + calibration_offset;
    DISTMASTER = compute_rssi_mismatch_maps( ...
        simulated_rssi_map, measured_rssi);

    output_file = fullfile( ...
        output_dir, sprintf("DISTMASTER27freq%d.mat", station_id));
    save(output_file, "DISTMASTER", "frequency_hz", ...
        "location_ids", "station_name", "-v7.3");

    fprintf("Saved %s: %d locations, %d frequencies.\n", ...
        output_file, num_locations, num_frequencies);
end
