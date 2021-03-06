% The classic function to calculate force due to fiber elasticity from a
% discrete collection of points, but this time with reference points that
% come from re-meshing. 
% Inputs: r = reference radius, (x,y) = current points on the contour, 
% gam = surface tension, k = tension constant, x0 = reference
% configurations
% Output: tension forces according to Section S2.1.1 in the supplementary
% info
function eforce = calcelasticforcerefs(r,x,y,gam,k,x0)
    N=length(x);
    eforce=zeros(N,2);
    h=2*pi*r/N;
    for iPt=1:N
        indexm1=iPt-1;
        if (indexm1==0)
            indexm1=N;
        end 
        indexp1=iPt+1;
        if (indexp1==N+1)
            indexp1=1;
        end 
        tao1=[x(iPt)-x(indexm1) y(iPt)-y(indexm1)]/norm(x0(iPt,:)-x0(indexm1,:));
        T1=gam+k*(norm(tao1)-1);
        tao1=tao1/norm(tao1);
        tao2=[x(indexp1)-x(iPt) y(indexp1)-y(iPt)]/norm(x0(iPt,:)-x0(indexp1,:));
        T2=gam+k*(norm(tao2)-1);
        tao2=tao2/norm(tao2);
        eforce(iPt,:)=(T2*tao2-T1*tao1)/h;
    end 
end 
        
    