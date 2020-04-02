% This function calculates the force on every ECM node. 
function fibforce=calcfibforce(curlocs,kteth,connections,reflocs)
    Nfib=length(connections);
    Npts=floor(length(curlocs)/Nfib+1e-2);
    % First is "pinning down force". This is handled by subtracting
    % all the tethered locations. 
    fibforce=curlocs-reflocs; % Subtract the reference locations. 
    for iFib=1:Nfib
        conlist=connections{iFib};
        for iCon=1:length(conlist)
            % Subtract the location of each tethered fiber (and its
            % associated place in the list). 
            jFib=conlist(iCon);
            fibforce((iFib-1)*Npts+1:iFib*Npts,:)=fibforce((iFib-1)*Npts+1:iFib*Npts,:)+...
                curlocs((iFib-1)*Npts+1:iFib*Npts,:)-curlocs((jFib-1)*Npts+1:jFib*Npts,:);
        end
    end
    % Multiply the result by -kteth. 
    fibforce=-kteth*fibforce; 
end    
    