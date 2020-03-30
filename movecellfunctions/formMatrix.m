% Function to compute the hydrodynamic matrix U = M*f using the method of
% regularized Stokeslets. Inputs: xks = points, eps = regularization
% parameter, mu = fluid viscositys
function M=formMatrix(xks,eps,mu)
    [nXk,~]=size(xks);
    M=zeros(2*nXk);
    for iE=1:nXk
        x=xks(iE,:);
        rk=(xks-x);
        nrk=sqrt(sum(rk.*rk,2));
        c = log(sqrt(nrk.^2+eps^2)+eps)-...
            eps*(sqrt(nrk.^2+eps^2)+2*eps)./((sqrt(nrk.^2+eps^2)+eps).*sqrt(nrk.^2+eps^2));
        d = (sqrt(nrk.^2+eps^2)+2*eps)./...
            ((sqrt(nrk.^2+eps^2)+eps).^2.*sqrt(nrk.^2+eps^2));
        m11 = 1/(4*pi*mu)*(-c+d.*rk(:,1).^2);
        m12 = 1/(4*pi*mu)*d.*rk(:,1).*rk(:,2);
        m21 = m12;
        m22 = 1/(4*pi*mu)*(-c+d.*rk(:,2).^2);
        M(2*iE-1,2*(1:nXk)-1) = m11;
        M(2*iE-1,2*(1:nXk)) = m12;
        M(2*iE,2*(1:nXk)-1) = m21;
        M(2*iE,2*(1:nXk)) = m22;
    end
end 