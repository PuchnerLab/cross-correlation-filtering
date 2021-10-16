function [REdges,Nautocounts,radii,Ncounts] = paircorr(coords, area, max_distance)
    kd = KDTreeSearcher(coords, 'BucketSize', 1);
    [~, d] = rangesearch(kd, coords, max_distance);
    distances = horzcat(d{:});
    N = size(coords, 1);
    density = N^2 / area;
    edges = [.01*ceil(min(distances)/.010):.030:.010*ceil(max(distances)/.010)];
    [counts, BinEdges] = histcounts(distances, edges);
    counts(1) = counts(1) - length(find(distances == 0));
    
    BinEdges=transpose(BinEdges);
    LEdges=BinEdges(1:end-1);
    REdges=BinEdges(2:end);
    BinWidth = BinEdges(2)-BinEdges(1);
    deltaAprime=pi*((LEdges+BinWidth).^2-(LEdges).^2);
    %radii=LEdges+((BinWidth)/2);
    radii = LEdges;
    deltaA=pi*2*radii.*(BinWidth)+pi.*(BinWidth).^2;
    Nfactor=deltaA.*N;
    counts=transpose(counts);
    Ncounts=counts./(2.*Nfactor);
    % figure
    % plot(radii,Ncounts)
    Nautocounts=counts./deltaAprime;

%     N=(XM_length);
%     density=(N.^2)/(40.96.^2);
%     edges = [.01*ceil(min(d_inside_ind2)/.010):.030:.010*ceil(max(d_inside_ind2)/.010)];
%     %edges = [0:.01:(.01*200)];
%     [counts,BinEdges]=histcounts(d_inside_ind2,edges);
%     counts(1) = counts(1) - length(find(d_inside_ind2 == 0));
%     BinEdges=transpose(BinEdges);
%     LEdges=BinEdges(1:end-1);
%     REdges=BinEdges(2:end);
%     BinWidth = BinEdges(2)-BinEdges(1);
%     deltaAprime=pi*((LEdges+BinWidth).^2-(LEdges).^2);
%     %radii=LEdges+((BinWidth)/2);
%     radii = LEdges;
%     deltaA=pi*2*radii.*(BinWidth)+pi.*(BinWidth).^2;
%     Nfactor=deltaA.*N;
%     counts=transpose(counts);
%     Ncounts=counts./(2.*Nfactor);
%     % figure
%     % plot(radii,Ncounts)
%     Nautocounts=counts./deltaAprime;
end