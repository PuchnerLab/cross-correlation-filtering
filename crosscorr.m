function [normcounts, centers] = crosscorr(coords1, coords2, max_distance, total_area)
    kd = KDTreeSearcher(coords1, 'BucketSize', 1);
    [~, d] = rangesearch(kd, coords2, max_distance);
    distances = horzcat(d{:});
    edges = 0:0.02:max_distance;
    [counts, BinEdges] = histcounts(distances, edges);
    areas = pi * diff(edges.^2);
    normcounts = total_area * counts ./ areas / size(coords1, 1) / size(coords2, 1);
    centers = (edges(1:end-1) + edges(2:end)) / 2;

%     density = (size(coords1, 1) * size(coords2, 1)) / total_area;
    
%     BinEdges=transpose(BinEdges);
%     LEdges=BinEdges(1:end-1);
%     REdges=BinEdges(2:end);
%     BinWidth = BinEdges(2)-BinEdges(1);
%     deltaAprime=pi*((LEdges+BinWidth).^2-(LEdges).^2);
%     radii=LEdges+((BinWidth)/2);
% 
%     radii_0=LEdges;
%     radii = radii_0;
% 
%     deltaA=pi*2*radii_0.*(BinWidth)+pi.*(BinWidth).^2;
% 
%     %deltaA=pi*2*radii.*(BinWidth);
%     Nfactor=deltaA.*density;
%     counts=transpose(counts);
%     Ncounts=counts./Nfactor;
end
