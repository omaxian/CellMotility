% Calculate the bending force in a contour with preferred curvature of a
% sphere. Inputs: (x,y) = contour locations, r = contour radius Kb =
% bending constant.
% Outputs: the bending force density
function bforce = calcbendingforce(x,y,r,Kb)
    N=length(x);
    h=2*pi*r/N;
    kap = 1/r;
    pts=[x y];
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
        Xss = 1/h^2*(pts(indexp1,:)+pts(indexm1,:)-2*pts(m,:));
        val=Kb/h^2*Xss*(1-kap/norm(Xss));
        bforce(indexp1,:)=bforce(indexp1,:)-val;
        bforce(m,:)=bforce(m,:)+2*val;
        bforce(indexm1,:)=bforce(indexm1,:)-val;
    end
end