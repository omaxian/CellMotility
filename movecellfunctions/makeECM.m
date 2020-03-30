% Make the ECM matrix. Outputs the centers of the fibers and the list of
% connections. 
% Inputs: Nfib = number of fibers, rcell = cell radius, nodebound =
% bounding box for the ECM nodes, spacing = min spacing between ECM points,
% eps = Stokeslet regularization parameter. 
function [fiberlocs, connectionsList] = makeECM(Nfib,rcell,nodebound,eps,spacing)
    fiberlocs=(rand(1,2)-[0 0.5])*nodebound; % first locations   
    rfib = 0;
    while (norm(fiberlocs(1,:)) < (rcell+rfib+2*eps))
        fiberlocs=(rand(1,2)-[0 0.5])*nodebound;
    end
    numfib=1;
    % Loop to make the rest of the fibers
    while (length(fiberlocs) < Nfib)
        proposed=(rand(1,2)-[0 0.5])*nodebound;
        diffs=proposed-fiberlocs;
        distances=(diffs(:,1).^2-diffs(:,2).^2).^(0.5);
        % The condition is that (1) the fibers are far enough from each other and 
        % (2) the fibers are initially far away from the cell so it can move 
        if (min(distances) > spacing && norm(proposed) > (rcell+rfib+2*eps))
            fiberlocs(numfib+1,:)=proposed; % add to the list if accepted
            numfib=numfib+1;
        end
    end
    % Do the neighbors list
    DT = delaunayTriangulation(fiberlocs);
    connectionslist={};
    for iFib=1:Nfib
        fibcons=[];
        for jFib=1:Nfib
            if (iFib ~=jFib)
                connect=isConnected(DT,iFib,jFib);
                if (connect)
                    fibcons=[fibcons jFib];
                end
            end
        end
        connectionsList{iFib}=fibcons;
    end
    triplot(DT);
end
    
            