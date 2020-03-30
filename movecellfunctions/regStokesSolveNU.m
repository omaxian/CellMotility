function [surfaceU,ns]=regStokesSolveNU(x,y,f,mu1,mu2,h,eps)
    Nm=length(x);
    lambda=mu2/mu1;
    kappa=(1-lambda)/(1+lambda);
    pvals=zeros(Nm,1);
    bterm=zeros(Nm,2);
    b=zeros(2*Nm,1);
    K=zeros(2*Nm);
    % Calculate the normal vectors
    ns=zeros(Nm,2);
    for iPt=1:Nm
        indexm1=iPt-1;
        if (indexm1==0)
            indexm1=Nm;
        end 
        indexp1=iPt+1;
        if (indexp1==Nm+1)
            indexp1=1;
        end 
        tao1=[x(iPt)-x(indexm1) y(iPt)-y(indexm1)]/h;
        tao2=[x(indexp1)-x(iPt) y(indexp1)-y(iPt)]/h;
        ns(iPt,:)=(tao2-tao1)/h;
        ns(iPt,:)=ns(iPt,:)/norm(ns(iPt,:));
    end 
    nxs=ns(:,1);
    nys=ns(:,2);
    for iPt=1:Nm
        xv=[x(iPt) y(iPt)];
        for jPt=1:Nm
            xk=[x(jPt) y(jPt)];
            r=xv-xk;
            rk=norm(r);
            fk=f(jPt,:);
            bterm(iPt,:)=bterm(iPt,:)-fk*(log(sqrt(rk^2+eps^2)+eps)-...
                eps*(sqrt(rk^2+eps^2)+2*eps)/((sqrt(rk^2+eps^2)+eps)*sqrt(rk^2+eps^2)))+...
                dot(fk,xv-xk)*(xv-xk)*((sqrt(rk^2+eps^2)+2*eps)/...
                ((sqrt(rk^2+eps^2)+eps)^2*sqrt(rk^2+eps^2)));
            if (norm(r)==0)
                K(2*iPt-1,2*jPt-1)=0;
                K(2*iPt-1,2*jPt)=0;
                K(2*iPt,2*jPt-1)=0;
                K(2*iPt,2*jPt)=0;
            else 
                K(2*iPt-1,2*jPt-1)=-4*r(1)*r(1)/norm(r)^4*(r(1)*nxs(jPt)+r(2)*nys(jPt));
                K(2*iPt-1,2*jPt)=-4*r(2)*r(1)/norm(r)^4*(r(1)*nxs(jPt)+r(2)*nys(jPt));
                K(2*iPt,2*jPt-1)=-4*r(1)*r(2)/norm(r)^4*(r(1)*nxs(jPt)+r(2)*nys(jPt));
                K(2*iPt,2*jPt)=-4*r(2)*r(2)/norm(r)^4*(r(1)*nxs(jPt)+r(2)*nys(jPt));
            end
        end
        bx0=1/(2*pi*mu1*(lambda+1))*bterm(iPt,:);
        b(2*iPt-1)=bx0(1);
        b(2*iPt)=bx0(2);
    end
   for iPt=1:Nm
        K(2*iPt-1,2*iPt-1)=-sum(K(2*iPt-1,1:2:2*Nm-1))-2*pi;
        K(2*iPt-1,2*iPt)=-sum(K(2*iPt-1,2:2:2*Nm));
        K(2*iPt,2*iPt-1)=-sum(K(2*iPt,1:2:2*Nm-1));
        K(2*iPt,2*iPt)=-sum(K(2*iPt,2:2:2*Nm))-2*pi;
   end
    K = -kappa*h/(2*pi)*K;
    W = (eye(2*Nm)-K)\b;
    % Rearrange
    surfaceU=zeros(Nm,2);
    for iPt=1:Nm
        surfaceU(iPt,1)=W(2*iPt-1);
        surfaceU(iPt,2)=W(2*iPt);
    end 
end 
            