% Function to select 2 protrusions on the surface of the cell for mechanism
% 2. Inputs: Nc = number of cortex points, points = if there are a
% pre-selected points for the protrusions
% Output: the 2 protrusions prot1 and prot2
% There are a lot of hacks in here to get the code to run. A good future
% project would be to make it more robust
function [prot1, prot2] = pickprots(Nc,points,t0,km)
    % Select the point for the protrusion
    % Prot 1 is ALWAYS BEHIND prot2, so that the loose part is between
    % prots 1 and 2. 
    prot1 = pickpoint(Nc,points(1));
    if (points(2) > -1)
        prot2 = points(2);
        return;
    end
    prot2=Nc/2;
    while (prot2 > Nc/4 && prot2 < 3*Nc/4)
        sgn=rand > 0.5;
        if (sgn==0)
            sgn=-1;
        end
        % Select prot2 approximately 30+-15 degrees from prot1
        % This is where the hacks are. For the first protrusion,we allow a
        % wider range of points. After that, the points have to be close
        % together so that the cell doesn't try to pass through ECM nodes. 
        if (t0)
            prot2 = prot1+sgn*(4+floor(rand*Nc/10));
            if ((prot1/Nc < 0.17 && prot2/Nc > 0.23) || ...
                (prot1/Nc > 0.23 && prot2/Nc < 0.17))
                prot2 = Nc/2;
            end
            if (km > 500) % if nucleus is stiff, space protrusions wider
               prot2 = prot1+sgn*(10+floor(rand*Nc/20));
            end
        else
            prot2 = prot1+sgn*(4+floor(rand*Nc/24));
            if (km > 500)
               prot2 = prot1+sgn*(10+floor(rand*Nc/20));
            end
        end
    end
    % Change which protusion is which and do proper modding
    if (prot2 < prot1) 
        temp=prot1;
        prot1=prot2;
        prot2=temp;
    end
    prot2 = mod(prot2,Nc); 
    prot2 = prot2 + (prot2==0)*Nc;
    prot1 = mod(prot1,Nc); 
    prot1 = prot1 + (prot1==0)*Nc;
end 