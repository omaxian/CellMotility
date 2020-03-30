function [prot1, prot2] = pickprots(Nc,point,numfibers,Npfib,ecmlocs,xc,yc)
    % Select the point for the protrusion
    % Prot 1 is ALWAYS BEHIND prot2, so that the loose part is between
    % prots 1 and 2. 
    prot1 = pickpoint(Nc,point,numfibers,Npfib,ecmlocs,xc,yc);
    prot2=Nc/2;
    while (prot2 > Nc/4 && prot2 < 3*Nc/4)
        sgn=rand > 0.5;
        if (sgn==0)
            sgn=-1;
        end
        prot2 = prot1+sgn*floor((pi/6+rand*pi/12)/(2*pi/Nc));
    end
    if (prot2 < prot1)
        temp=prot1;
        prot1=prot2;
        prot2=temp;
    end
    if (prot2 > Nc)
        prot2=prot2-Nc;
    elseif (prot2 < 1)
        prot2=prot2+Nc;
    end
    if (prot1 > Nc)
        prot1=prot1-Nc;
    elseif (prot1 < 1)
        prot1=prot1+Nc;
    end
    
end 