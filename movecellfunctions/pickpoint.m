% This function selects the point where the protrusion will form.
function jPt = pickpoint(N,jPtin)
    if (jPtin < 0)
        jPt=floor((rand-0.5)*N/2)+1; % Bias protrusions in the front. 
        if (jPt < 1) % mod by N
            jPt=jPt+N;
        end
    else 
        jPt=jPtin; %there is also an option for the user to just pass in the point. 
    end
end
