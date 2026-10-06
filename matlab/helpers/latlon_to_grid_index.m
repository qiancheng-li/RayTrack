function grid_indices = latlon_to_grid_index( ...
    latlon, longitude_step, latitude_step, min_longitude, min_latitude)
%LATLON_TO_GRID_INDEX Convert latitude/longitude coordinates to grid indices.

    if nargin < 5
        longitude_step = 4.008212542316869e-06;
        latitude_step = 2.699970153230424e-06;
        min_longitude = -122.3105221345138;
        min_latitude = 47.651033679860470;
    end

    grid_indices = (latlon - [min_latitude, min_longitude]) ./ ...
        [latitude_step, longitude_step] + 1;
end
