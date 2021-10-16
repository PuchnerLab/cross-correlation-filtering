function [cluster_filter, coord_idx, coords1_coloc] = separate_clusters(mlist1, mlist2, distance, minimumNumber1, minimumNumber2)
    % Check if a cluster colocalizes with another cluster. Returns the
    % indices of mlist1.blinking.ClusterList that mlist2 colocalizes
    % with.
    mlist1_idx = 1:length(mlist1.blinking.newx);
    kd = KDTreeSearcher([mlist2.blinking.newx(:), mlist2.blinking.newy(:)]);
    [idx, dist] = rangesearch(kd, [mlist1.blinking.newx(:), mlist1.blinking.newy(:)], distance);
    coords2_idx = cellfun(@(d, i) i((d <= distance) & (length(i) >= minimumNumber2)), dist, idx, 'UniformOutput', false);
    % coords2_idx = cellfun(@(d, i) i((d <= distance)), dist, idx, 'UniformOutput', false);
    coords1_coloc = ~cellfun(@isempty, coords2_idx);
    cluster1_select = cellfun(@(y) sum(ismember(y, mlist1_idx(coords1_coloc))) > 0,...
        mlist1.blinking.ClusterList);
    cluster_filter = cluster1_select & (mlist1.blinking.clusternum >= minimumNumber1);
    % cluster_filter = cluster1_select;
    coord_idx = sort(cell2mat(mlist1.blinking.ClusterList(cluster_filter)));
end
