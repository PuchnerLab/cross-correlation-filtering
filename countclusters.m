function MList = countclusters(MList, countradius, excludesingleblinkers)

%% Check if the structure contains the molecule list under a .master field.
% If it doesn't exist an error message will appear notifying useer.
    if isfield(MList,'blinking') == 0 
       warning('Structure field "blinking" containing a blink corrected molecules list does not exist in "%s",\n Check the name or run "%s = ExtractBlinking" first.',...
           inputname(1),inputname(1));
       return;
    end

%% Get required files from MList.
    
%   The new components that will be added to the molecules list
    countx = MList.blinking.newx;
    county = MList.blinking.newy;
    Lifetime = MList.blinking.Lifetime;
    
   % clusterid = zeros(length(countx),1);
    clusternum = ones(1,length(countx));
    
    numMol = length(countx);
    index = 1:numMol;
    ClusterList = cell(1,numMol);
        for k = 1:numMol
            ClusterList{k} = index(k);
        end
    
    if(excludesingleblinkers == 1)
        exclusion = find(Lifetime == 1);
        countx(exclusion) = [];
        county(exclusion) = [];
        Lifetime(exclusion) = [];
        clusternum(exclusion) = [];
        ClusterList(exclusion) = [];
    end
        
    n=1;
    while n<length(countx)
        
        for k = 1:length(countx)
            if k~=n
            if sqrt((countx(k)-countx(n))^2 + (county(k)-county(n))^2)<countradius
                %Identifying molecules that can belong to a cluster.
                countx(n) = (clusternum(n)*countx(n) + clusternum(k)*countx(k))/(clusternum(n)+clusternum(k));     
                county(n) = (clusternum(n)*county(n) + clusternum(k)*county(k))/(clusternum(n)+clusternum(k));
                clusternum(n) = clusternum(n)+clusternum(k);
                
                ClusterList{n} = [ClusterList{n} index(k)];
                ClusterList(k) = [];
                
                %Deleting stuff.
                countx(k) = [];
                county(k) = [];
                clusternum(k) = [];
                index(k) = [];
                
                n=n-1; %go back one to restart calculation from that new weighted location.
                break  %Exit for statement
                
            end
            end
        end
        n=n+1;
    end
    
    MList.blinking.countx = countx;
    MList.blinking.county = county;
    MList.blinking.clusternum = clusternum;
    MList.blinking.ClusterList = ClusterList;
                

    
end