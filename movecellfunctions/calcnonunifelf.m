% Calculate the elastic force when there is a reference confgiguration
% (like after the cell has been re-meshed) AND a non-uniform tension
% constant.
% Inputs: (x,y) = current points on the contour, r = reference radius, 
% Kb = bending constant in the regions where it exists (elsewhere it is 0), 
% prot1 = where the first protrusion is bound, prot2 = where the second
% protrusion is bound, (x0,y0) = reference configuration. 
% (gamsoft, ksoft) = tension constants for the soft part of the cortex
% (gamhard,khard)= tension constants for the hard (contracting) part of the cortex
% We assign K_b to be Kb in the region behind prot1 and prot2 (inclusive), 
% otherwise it's zero.
% Output: elastic forces according to Section S2.1.1 in the supplementary
% info
function eforce = calcnonunifelf(r,x,y,x0,y0,gamsoft,gamhard,ksoft,khard,prot1,prot2)
    nrl=0.1; % the rest length is set to 0.1 in the hard region
    N=length(x);
    eforce=zeros(N,2);
    h=2*pi*r/N;
    if (prot1-prot2 > 0) % opposite sides of the 0 marker
        looseinds=[prot1+1:N 1:prot2-1];
    else
        looseinds=prot1+1:prot2-1;
    end
    % we need soft in between because they're on the same side of the 1st point
    for iPt=1:N
        indexm1=iPt-1;
        if (indexm1==0)
            indexm1=N;
        end 
        indexp1=iPt+1;
        if (indexp1==N+1)
            indexp1=1;
        end 
        if (iPt~=prot1 && iPt~=prot2)
            if (sum(looseinds==iPt) > 0) % iPt is in the soft section
                tao1=[x(iPt)-x(indexm1) y(iPt)-y(indexm1)]/...
                    norm([x0(iPt)-x0(indexm1) y0(iPt)-y0(indexm1)]);
                T1=gamsoft+ksoft*(norm(tao1)-1);
                tao1=tao1/norm(tao1);
                tao2=[x(indexp1)-x(iPt) y(indexp1)-y(iPt)]/...
                    norm([x0(iPt)-x0(indexp1) y0(iPt)-y0(indexp1)]);
                T2=gamsoft+ksoft*(norm(tao2)-1);
                tao2=tao2/norm(tao2);
                eforce(iPt,:)=(T2*tao2-T1*tao1)/h;
            else % iPt is in the hard section
                tao1=[x(iPt)-x(indexm1) y(iPt)-y(indexm1)]/...
                    norm([x0(iPt)-x0(indexm1) y0(iPt)-y0(indexm1)]);
                T1=gamhard+khard*(norm(tao1)-nrl);
                tao1=tao1/norm(tao1);
                tao2=[x(indexp1)-x(iPt) y(indexp1)-y(iPt)]/...
                    norm([x0(iPt)-x0(indexp1) y0(iPt)-y0(indexp1)]);
                T2=gamhard+khard*(norm(tao2)-nrl);
                tao2=tao2/norm(tao2);
                eforce(iPt,:)=(T2*tao2-T1*tao1)/h;
            end
        elseif (iPt==prot1) % it is prot1
            tao1=[x(iPt)-x(indexm1) y(iPt)-y(indexm1)]/...
                    norm([x0(iPt)-x0(indexm1) y0(iPt)-y0(indexm1)]);
            T1=gamhard+khard*(norm(tao1)-nrl); % hard behind it
            tao1=tao1/norm(tao1);
            tao2=[x(indexp1)-x(iPt) y(indexp1)-y(iPt)]/...
                    norm([x0(iPt)-x0(indexp1) y0(iPt)-y0(indexp1)]);
            T2=gamsoft+ksoft*(norm(tao2)-1); % soft ahead of it
            tao2=tao2/norm(tao2);
            eforce(iPt,:)=(T2*tao2-T1*tao1)/h;
        else % it is prot2
            tao1=[x(iPt)-x(indexm1) y(iPt)-y(indexm1)]/...
                    norm([x0(iPt)-x0(indexm1) y0(iPt)-y0(indexm1)]);
            T1=gamsoft+ksoft*(norm(tao1)-1); % soft behind it
            tao1=tao1/norm(tao1);
            tao2=[x(indexp1)-x(iPt) y(indexp1)-y(iPt)]/...
                    norm([x0(iPt)-x0(indexp1) y0(iPt)-y0(indexp1)]);
            T2=gamhard+khard*(norm(tao2)-nrl); % hard ahead of it
            tao2=tao2/norm(tao2);
            eforce(iPt,:)=(T2*tao2-T1*tao1)/h;
        end
    end
end 
        
    