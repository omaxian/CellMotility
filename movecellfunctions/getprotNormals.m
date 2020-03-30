% Compute the outward pointing normal vector at a point on a contour. 
% Inputs: prot1 = point where we want the normal, (xc,yc) = points on the
% contour
% Outputs: the normal vector at that points
function np1=getprotNormals(prot1,xc,yc)
    Nc=length(xc);
    % Normal vectors
    iPt=prot1;
    indexm1=iPt-1;
    if (indexm1==0)
        indexm1=Nc;
    end 
    indexp1=iPt+1;
    if (indexp1==Nc+1)
        indexp1=1;
    end 
    tao1=[xc(iPt)-xc(indexm1) yc(iPt)-yc(indexm1)];
    tao1=tao1/norm(tao1);
    tao2=[xc(indexp1)-xc(iPt) yc(indexp1)-yc(iPt)];
    tao2=tao2/norm(tao2);
    np1=(tao2-tao1);
    np1=np1/norm(np1);
    if (dot([xc(prot1) yc(prot1)]-mean([xc yc]),np1) < 0)
        np1=-np1;
    end
end