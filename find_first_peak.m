function [loc_percentile, idx] = find_first_peak(x, minwidth)
    [~, locs] = findpeaks(-x, 'minpeakwidth', minwidth);
    loc = locs(1);
    [~, idx] = min(abs(x(1:loc) - (x(1) + x(loc)) / 2));
    percentile = 0.95;
    [~, loc_percentile] = min(abs(x(1:loc) - (x(1) + x(loc)) * (1 - percentile)));
end
