function [cluster_filters, coord_idx, coords_coloc, included, excluded] = pipeline(coords1, coords2, fovArea, maxDistance, countclusters_distance, minimumNumber, units)
    % % Calculate pair correlation to determine cluster assignment distance.

%     % Parameter
%     fovArea = 40.96^2;

%     minimumNumber = 2;  % minimum number of localizations needed to classify as a cluster.

%     % Calculate automatically
%     countclusters_distance = 0.4;  % um

%     maxDistance = 1; % um. Maximum distance in pair-correlation and cross-correlation calculations.

    % % Determine cluster separation
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

    % % Angel's code; Calculate cross-correlation to determine co-localization cutoff distance.
    [normcounts, binCenters] = crosscorr(...
        coords1,...
        coords2,...
        maxDistance,...
        fovArea);

%     [~, optimal_index] = min(abs(normcounts - ((normcounts(1) - normcounts(end)) / 2 + normcounts(end))));
%     optimal_distance = binCenters(optimal_index); % bigr(argmax_both); % bigr(max(argmax_pop1, argmax_pop2));

    [~, half_index] = min(abs(normcounts - ((normcounts(1) - normcounts(end)) / 2 + normcounts(end))));
    half_distance = binCenters(half_index); % bigr(argmax_both); % bigr(max(argmax_pop1, argmax_pop2));

    [~, bottom_index] = min(abs(normcounts - ((normcounts(1) - normcounts(end)) * 0.01 + normcounts(end))));
    bottom_distance = binCenters(bottom_index); % bigr(argmax_both); % bigr(max(argmax_pop1, argmax_pop2));

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


    % % Assign clusters
    coords1_mlist = add_blinking_fields(coords1);
    coords1_mlist = countclusters(coords1_mlist, countclusters_distance, 0);
    coords1_mlist = compute_cluster_size(coords1_mlist);

    coords2_mlist = add_blinking_fields(coords2);
    coords2_mlist = countclusters(coords2_mlist, countclusters_distance, 0);
    coords2_mlist = compute_cluster_size(coords2_mlist);

%     % % Dushyant's code; Calculate cross-correlation to determine co-localization cutoff distance.
%     [REdgesbig, Nautocountsbig, radiibig, Ncountsbig, xym1in, xym2in, bigr]...
%         = memorycrosscorrelationwithmoleculeseparation_direct_data(...
%         [coords1_mlist.blinking.newx, coords1_mlist.blinking.newy], ...
%         [coords2_mlist.blinking.newx, coords2_mlist.blinking.newy],...
%         maxDistance);

    % Filter colocalized clusters based on stoichiometry. The filtering is
    % done both ways to support stoichiometry thresholds for each of the
    % two populations.
    cluster_filters = cell([1 2]);
    coord_idx = cell([1 2]);
    coords_coloc = cell([1 2]);
    [cluster_filters{1}, coord_idx{1}, coords_coloc{1}] = separate_clusters(...
        coords1_mlist, coords2_mlist, optimal_distance,...
        minimumNumber, minimumNumber);
    [cluster_filters{2}, coord_idx{2}, coords_coloc{2}] = separate_clusters(...
        coords2_mlist, coords1_mlist, optimal_distance,...
        minimumNumber, minimumNumber);

    included = cell([1 2]);
    included{1} = coords1(coords_coloc{1}, :);
    included{2} = coords2(coords_coloc{2}, :);

    excluded = cell([1 2]);
    excluded{1} = coords1(~coords_coloc{1}, :);
    excluded{2} = coords2(~coords_coloc{2}, :);

%     included = cell([1 2]);
%     included{1} = [coords1_mlist.blinking.newx(coord_idx{1}),...
%                    coords1_mlist.blinking.newy(coord_idx{1})];
%     included{2} = [coords2_mlist.blinking.newx(coord_idx{2}),...
%                    coords2_mlist.blinking.newy(coord_idx{2})];
% %     included1 = [coords1_mlist.blinking.newx(coord_idx{1}),...
% %                  coords1_mlist.blinking.newy(coord_idx{1})];
% %     included2 = [coords2_mlist.blinking.newx(coord_idx{2}),...
% %                  coords2_mlist.blinking.newy(coord_idx{2})];
% 
%     excluded = cell([1 2]);
%     excluded{1} = coords1(~ismember(coords1, included{1}, 'rows'), :);
%     excluded{2} = coords2(~ismember(coords2, included{2}, 'rows'), :);
%     % excluded1 = coords1(~ismember(coords1, included1, 'rows'), :);
%     % excluded2 = coords2(~ismember(coords2, included2, 'rows'), :);
end
