% Calculate the bending force when there is a reference confgiguration
% (like after the cell has been re-meshed) AND a non-uniform bending
% constant.
% Inputs: (x,y) = current points on the contour, r = reference radius, 
% Kb = bending constant in the regions where it exists (elsewhere it is 0), 
% prot1 = where the first protrusion is bound, prot2 = where the second
% protrusion is bound, x0 = reference configuration. 
% We assign K_b to be Kb in the region behind prot1 and prot2 (inclusive), 
% otherwise it's zero.
% Output: bending forces according to Section S2.1.2 in the supplementary
% info
function bforce = calcnonunifbend(x,y,r,x0,Kb,prot1,prot2)
    N=length(x);
    kap = 1/r;
    pts=[x y];
    h=2*pi*r/N;
    if (prot1-prot2 > 0) % opposite sides of the 0 marker
        looseinds=[prot1+1:N 1:prot2-1];
    else
        looseinds=prot1+1:prot2-1;
    end
    bforce=zeros(N,2);
    for m=1:N
        indexm1=m-1;
        if (indexm1==0)
            indexm1=N;
        end 
        indexp1=m+1;
        if (indexp1==N+1)
            indexp1=1;
        end 
        if (m~=prot1 && m~=prot2 && sum(looseinds==m) > 0)
            % Soft section, do nothing
        else
            Xss = 2/norm(x0(indexp1,:)-x0(indexm1,:))*...
            ((pts(indexp1,:)-pts(m,:))/norm(x0(indexp1,:)-x0(m,:))+...
            (-pts(m,:)+pts(indexm1,:))/norm(x0(indexm1,:)-x0(m,:)));
            val=Kb/h*Xss*(1-kap/norm(Xss));
            bforce(indexp1,:)=bforce(indexp1,:)-val/norm(x0(indexp1,:)-x0(m,:));
            bforce(m,:)=bforce(m,:)+val/norm(x0(indexp1,:)-x0(m,:))+val/norm(x0(indexm1,:)-x0(m,:));
            bforce(indexm1,:)=bforce(indexm1,:)-val/norm(x0(indexm1,:)-x0(m,:));
        end
    end
end 
        
    