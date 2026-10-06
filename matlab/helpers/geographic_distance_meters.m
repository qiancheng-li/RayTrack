function distance_meters = geographic_distance_meters( ...
    latitude, longitude, reference_latitude, reference_longitude)
%GEOGRAPHIC_DISTANCE_METERS Approximate point-to-point geographic distance.

    [x_meters, y_meters] = latlon_to_local_xy( ...
        latitude, longitude, reference_latitude, reference_longitude);
    distance_meters = hypot(x_meters, y_meters);
end
