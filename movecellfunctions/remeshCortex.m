function [x_new, x_new_ref,newprot1,newprot2,np1,np2] = remeshCortex(Nc,x_now,x0,prot1,prot2)
    N=Nc;
    x_new=x_now;
    x_new_ref=x0;
    % The old interpolant
%     shifted_0=[x0(1,:); x0(1:end-1,:)];
%     d_0=sqrt((x0(:,1)-shifted_0(:,1)).^2+(x0(:,2)-shifted_0(:,2)).^2);
%     d_cum0=d_0;
%     for iP=1:N
%         d_cum0(iP)=sum(d_0(1:iP));
%     end
    % The new interpolant
    shifted=[x_now(1,:); x_now(1:end-1,:)];
    d=sqrt((x_now(:,1)-shifted(:,1)).^2+(x_now(:,2)-shifted(:,2)).^2);
    d_cum=d;
    for iP=1:N
        d_cum(iP)=sum(d(1:iP));
    end
    dlast=d_cum(end)+norm(x_now(end,:)-x_now(1,:));
    % Sample the new at equal distances
    hnew=(dlast-d_cum(1))/N;
    dsamps=[d_cum(1):hnew:dlast-hnew]';
    [~,replace1]=min(abs(dsamps-d_cum(prot1)));
    dsamps(replace1)=d_cum(prot1);
    [~,replace2]=min(abs(dsamps-d_cum(prot2)));
    dsamps(replace2)=d_cum(prot2);
    % Get the s values from the d values
    for iP=1:N
        if (dsamps(iP)>d_cum(end))
            d_cum(1)=dlast;
        else
            d_cum(1)=0;
        end
        negds=d_cum-dsamps(iP);
        negds(negds>1e-10)=1000;
        [away1,ind1]=min(abs(negds));
        temp=d_cum(ind1);
        d_cum(ind1)=1000;
        posds=d_cum-dsamps(iP);
        posds(posds<0)=1000;
        [away2,ind2]=min(abs(posds));
        d_cum(ind1)=temp;
        node1=min(ind1,ind2);
        node2=max(ind1,ind2);
        dbtw=d_cum(node2)-d_cum(node1);
        if (dbtw < 0)
            node1=max(ind1,ind2);
            node2=min(ind1,ind2);
            dbtw=-dbtw;
        end
        if (node1==ind1)
            d1=away1;
        else
            d1=away2;
        end
        btw=d1/dbtw;
        if (isnan(btw))
            btw=0;
        end
        % Get the corresponding pt
        x_new(iP,:)=x_now(node1,:)+btw*(x_now(node2,:)-x_now(node1,:));
        nodenum(iP)=node1+btw;
        % Get the corresponding pt in the reference configuration!
        x_new_ref(iP,:)=x0(node1,:)+btw*(x0(node2,:)-x0(node1,:));
    end
    % Re-order so max x is in front
    [~,ind]=max(x_new(:,1));
    x_new=[x_new(ind:end,:); x_new(1:ind-1,:)];
    x_new_ref=[x_new_ref(ind:end,:); x_new_ref(1:ind-1,:)];
    % Extra info
    xc=x_now(:,1);
    yc=x_now(:,2);
    newxs=x_new(:,1);
    newys=x_new(:,2);
    [val,newprot1]=min((newxs-xc(prot1)).^2+(newys-yc(prot1)).^2);
    [val2,newprot2]=min((newxs-xc(prot2)).^2+(newys-yc(prot2)).^2);
    if (val > 1e-12 || val2 > 1e-12)
        keyboard
    end
    if (newprot1 == newprot2)
        keyboard
    end
    scatter(newxs,newys)
    % Normal vectors
    np1=getprotNormals(newprot1,newxs,newys);
    np2=getprotNormals(newprot2,newxs,newys);
%     iPt=newprot1;
%     indexm1=iPt-1;
%     if (indexm1==0)
%         indexm1=Nc;
%     end 
%     indexp1=iPt+1;
%     if (indexp1==Nc+1)
%         indexp1=1;
%     end 
%     tao1=[newxs(iPt)-newxs(indexm1) newys(iPt)-newys(indexm1)];
%     tao1=tao1/norm(tao1);
%     tao2=[newxs(indexp1)-newxs(iPt) newys(indexp1)-newys(iPt)];
%     tao2=tao2/norm(tao2);
%     np1=(tao2-tao1);
%     np1=np1/norm(np1);
%     if (dot([xc(prot1) yc(prot1)]-mean([xc yc]),np1) > 0)
%         np1=-np1;
%     end
%     iPt=newprot2;
%     indexm1=iPt-1;
%     if (indexm1==0)
%         indexm1=Nc;
%     end 
%     indexp1=iPt+1;
%     if (indexp1==Nc+1)
%         indexp1=1;
%     end 
%     tao1=[newxs(iPt)-newxs(indexm1) newys(iPt)-newys(indexm1)];
%     tao1=tao1/norm(tao1);
%     tao2=[newxs(indexp1)-newxs(iPt) newys(indexp1)-newys(iPt)];
%     tao2=tao2/norm(tao2);
%     np2=(tao2-tao1);
%     np2=np2/norm(np2);
%     if (dot([xc(prot2) yc(prot2)]-mean([xc yc]),np2) > 0)
%         np2=-np2;
%     end
end
