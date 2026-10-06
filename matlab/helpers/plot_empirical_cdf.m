function line_handle = plot_empirical_cdf(values)
%PLOT_EMPIRICAL_CDF Plot the empirical cumulative distribution function.

    values = values(isfinite(values));
    values = sort(values(:));

    if isempty(values)
        error("plot_empirical_cdf:NoFiniteValues", ...
            "The input must contain at least one finite value.");
    end

    cumulative_probability = (1:numel(values)).' / numel(values);
    plot_x = [values(1); values];
    plot_y = [0; cumulative_probability];

    line_handle = stairs(plot_x, plot_y, "LineWidth", 1.5);
end
