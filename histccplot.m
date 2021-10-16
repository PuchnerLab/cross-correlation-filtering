function [REdges,Nautocounts,radii,Ncounts] = histccplot(atg_length,XM_length,d_inside_ind2,smallr1,big_r1,matfilename)
tic
disp('plot started')
N=(atg_length*XM_length);
density=N/(40.96^2);

edges = [.05*ceil(min(d_inside_ind2)/.050):.050:.050*ceil(max(d_inside_ind2)/.050)];

[counts,BinEdges]=histcounts(d_inside_ind2,edges);
      
        BinEdges=transpose(BinEdges);
        LEdges=BinEdges(1:end-1);
        REdges=BinEdges(2:end);
        BinWidth = BinEdges(2)-BinEdges(1);
        deltaAprime=pi*((LEdges+BinWidth).^2-(LEdges).^2);
        radii=LEdges+((BinWidth)/2);
        
                
        radii_0=LEdges;
        radii = radii_0;

        deltaA=pi*2*radii_0.*(BinWidth)+pi.*(BinWidth).^2;
        
        
        
        %deltaA=pi*2*radii.*(BinWidth);
        Nfactor=deltaA.*density;
        counts=transpose(counts);
        Ncounts=counts./Nfactor;
% % figure
% % plot(radii,Ncounts)
Nautocounts=counts./deltaAprime;
titlestr = [matfilename ' ' num2str(smallr1), ' to ' ,num2str(big_r1), ' new.jpg'];
titlestr2 = [matfilename ' ' num2str(smallr1), ' to ' ,num2str(big_r1), ' new.fig'];

saveas(gcf,titlestr)
saveas(gcf,titlestr2)
%Nautocounts2 = gather(Nautocounts);
% % figure
% % plot(REdges,Nautocounts)
titlestr3 = [matfilename ' ' num2str(smallr1), ' to ' ,num2str(big_r1), ' new counts.jpg'];
titlestr4 = [matfilename ' ' num2str(smallr1), ' to ' ,num2str(big_r1), ' new counts.fig'];
saveas(gcf,titlestr3)
saveas(gcf,titlestr4)
toc
disp('plot ended')
filename2 = [matfilename num2str(smallr1), ' to ' ,num2str(big_r1) ' .mat'];
save(filename2,'REdges','Nautocounts','radii','Ncounts','-v7.3')

%clear BinEdges radii Ncounts counts Nautocounts deltaAprime BinWidth REdges LEdges BinEdges
end