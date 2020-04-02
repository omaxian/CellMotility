% Regularized Stokeset solve (SLOW using for loops). 
% Inputs: (xe,ye) = evaluation points, (xm,ym) = points where force is, 
% f = force in pN/um at those points, eps = regularization parameter, 
% mu = fluid viscosity
function [pvals,uvals]=regStokesSolve(xe,ye,xm,ym,f,eps,mu)
    Neval=length(xe);
    Nm=length(xm);
    pvals=zeros(Neval,1);
    uvals=zeros(Neval,2);
    for iPt=1:Neval
        x=[xe(iPt) ye(iPt)];
        for jPt=1:Nm
            xk=[xm(jPt) ym(jPt)];
            rk=norm(x-xk);
            fk=f(jPt,:);
%             pvals(iPt)=pvals(iPt)+1/(2*pi)*dot(fk,x-xk)*...
%                 (rk^2+2*eps^2)/(rk^2+eps^2)^2;
%             uvals(iPt,:)=uvals(iPt,:)+1/(4*pi*mu)*(-fk*(1/2*log(rk^2+eps^2) - ...
%                 eps^2/(rk^2+eps^2)) + dot(fk,x-xk)*(x-xk)*1/(rk^2+eps^2));
            pvals(iPt)=pvals(iPt)+1/(2*pi)*dot(fk,x-xk)*...
                ((rk^2+2*eps^2+eps*sqrt(rk^2+eps^2))/...
                ((sqrt(rk^2+eps^2)+eps)*(rk^2+eps^2)^1.5));
            uvals(iPt,:)=uvals(iPt,:)-fk/(4*pi*mu)*(log(sqrt(rk^2+eps^2)+eps)-...
                eps*(sqrt(rk^2+eps^2)+2*eps)/((sqrt(rk^2+eps^2)+eps)*sqrt(rk^2+eps^2)))+...
                1/(4*pi*mu)*dot(fk,x-xk)*(x-xk)*((sqrt(rk^2+eps^2)+2*eps)/...
                ((sqrt(rk^2+eps^2)+eps)^2*sqrt(rk^2+eps^2)));
        end
    end 
            