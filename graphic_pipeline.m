function graphic_pipeline(coords1, coords2, maxDistance, fovArea, units)
    % % Calculate pair correlation to determine cluster assignment distance.
    [~, ~, radii1, Ncounts1] = paircorr(coords1, fovArea, maxDistance);
    [~, ~, radii2, Ncounts2] = paircorr(coords2, fovArea, maxDistance);

    [~, pair_index1] = min(abs(Ncounts1 - ((Ncounts1(1) - Ncounts1(end)) * 0.01 + Ncounts1(end))));
    pair_distance1 = radii1(pair_index1);
    [~, pair_index2] = min(abs(Ncounts2 - ((Ncounts2(1) - Ncounts2(end)) * 0.01 + Ncounts2(end))));
    pair_distance2 = radii1(pair_index2);

    figure
    hold on
    plot(radii1, Ncounts1)
    plot([pair_distance1, pair_distance1], [min(Ncounts1), max(Ncounts1)]);
    text(pair_distance1, max(Ncounts1), num2str(pair_distance1))
    hold off
    xlabel(['Distance ' '(' units ')'])
    ylabel(['Radial distribution 1 ' '(#/' units '^2)'])
    
    figure
    hold on
    plot(radii2, Ncounts2)
    plot([pair_distance2, pair_distance2], [min(Ncounts2), max(Ncounts2)]);
    text(pair_distance2, max(Ncounts2), num2str(pair_distance2))
    hold off
    xlabel(['Distance ' '(' units ')'])
    ylabel(['Radial distribution 2 ' '(#/' units '^2)'])

    % % Calculate cross-correlation to determine co-localization cutoff distance.
    [normcounts, binCenters] = crosscorr(...
        coords1,...
        coords2,...
        maxDistance,...
        fovArea);

    [~, half_index] = min(abs(normcounts - ((normcounts(1) - normcounts(end)) / 2 + normcounts(end))));
    half_distance = binCenters(half_index);

    [~, bottom_index] = min(abs(normcounts - ((normcounts(1) - normcounts(end)) * 0.01 + normcounts(end))));
    bottom_distance = binCenters(bottom_index);

    figure
    hold on
    plot(binCenters, normcounts)
    plot([half_distance, half_distance], [min(normcounts), max(normcounts)]);
    text(half_distance, max(normcounts), num2str(half_distance))
    plot([bottom_distance, bottom_distance], [min(normcounts), max(normcounts)]);
    text(bottom_distance, max(normcounts), num2str(bottom_distance))
    hold off
    xlabel(['Distance ' '(' units ')'])
    ylabel('Cross-correlation')
end
