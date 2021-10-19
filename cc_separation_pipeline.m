function  [coords_coloc, included, excluded, cluster_filters, coord_idx] = cc_separation_pipeline(coords1, coords2, countclusters_distances, separation_distance, minimumNumbers)
    % % Assign clusters
    coords1_mlist = add_blinking_fields(coords1);
    coords1_mlist = countclusters(coords1_mlist, countclusters_distances(1), 0);
    coords1_mlist = compute_cluster_size(coords1_mlist);

    coords2_mlist = add_blinking_fields(coords2);
    coords2_mlist = countclusters(coords2_mlist, countclusters_distances(2), 0);
    coords2_mlist = compute_cluster_size(coords2_mlist);

    % Filter colocalized clusters based on stoichiometry. The filtering is
    % done both ways to support stoichiometry thresholds for each of the
    % two populations.
    cluster_filters = cell([1 2]);
    coord_idx = cell([1 2]);
    coords_coloc = cell([1 2]);
    [cluster_filters{1}, coord_idx{1}, coords_coloc{1}] = separate_clusters(...
        coords1_mlist, coords2_mlist, separation_distance,...
        minimumNumbers(1), minimumNumbers(2));
    [cluster_filters{2}, coord_idx{2}, coords_coloc{2}] = separate_clusters(...
        coords2_mlist, coords1_mlist, separation_distance,...
        minimumNumbers(2), minimumNumbers(1));

    included = cell([1 2]);
    included{1} = coords1(coords_coloc{1}, :);
    included{2} = coords2(coords_coloc{2}, :);

    excluded = cell([1 2]);
    excluded{1} = coords1(~coords_coloc{1}, :);
    excluded{2} = coords2(~coords_coloc{2}, :);
end
