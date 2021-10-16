function mlist = compute_cluster_size(mlist)
    mlist.blinking.clustersize = zeros(size(mlist.blinking.countx));
    for ii = 1:numel(mlist.blinking.countx)
        idx = mlist.blinking.ClusterList{ii};
        xs = mlist.blinking.newx(idx);
        ys = mlist.blinking.newy(idx);
        centerx = mlist.blinking.countx(ii);
        centery = mlist.blinking.county(ii);
        radii = sqrt((centerx - xs).^2 + (centery - ys).^2);
        mlist.blinking.clustersize(ii) = mean(radii);
    end
end
