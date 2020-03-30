function [x_new, x_new_ref] = remeshNucleus(Nc,x_now,x0)
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
end
