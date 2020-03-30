% The classic function to calculate force due to fiber elasticity from a
% discrete collection of points. 
% Inputs: r = reference config radius, (x,y) = contour points, gam =
% surface tension, k = spring constant
function eforce = calcelasticforce(r,x,y,gam,k)
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
        tao1=[x(iPt)-x(indexm1) y(iPt)-y(indexm1)]/h;
        T1=gam+k*(norm(tao1)-1);
        tao1=tao1/norm(tao1);
        tao2=[x(indexp1)-x(iPt) y(indexp1)-y(iPt)]/h;
        T2=gam+k*(norm(tao2)-1);
        tao2=tao2/norm(tao2);
        eforce(iPt,:)=(T2*tao2-T1*tao1)/h;
    end 
end 
        
    