function [x_meters, y_meters] = latlon_to_local_xy( ...
    latitude, longitude, reference_latitude, reference_longitude)
%LATLON_TO_LOCAL_XY Convert coordinates to local east/north distances.
%
% This equirectangular approximation is suitable for the small geographic
% area covered by the experiment.

    earth_radius_meters = 6371000;
    latitude_radians = deg2rad(latitude);
    longitude_radians = deg2rad(longitude);
    reference_latitude_radians = deg2rad(reference_latitude);
    reference_longitude_radians = deg2rad(reference_longitude);

    latitude_delta = latitude_radians - reference_latitude_radians;
    longitude_delta = longitude_radians - reference_longitude_radians;
    mean_latitude = (latitude_radians + reference_latitude_radians) / 2;

    x_meters = earth_radius_meters .* longitude_delta .* cos(mean_latitude);
    y_meters = earth_radius_meters .* latitude_delta;
end
