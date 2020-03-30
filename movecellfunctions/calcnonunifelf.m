function eforce = calcnonunifelf(r,x,y,x0,y0,gamsoft,gamhard,ksoft,khard,prot1,prot2)
    nrl=0.1;
    N=length(x);
    eforce=zeros(N,2);
    h=2*pi*r/N;
    if (prot1-prot2 > 0) % opposite sides of the 0 marker
        looseinds=[prot1+1:N 1:prot2-1];
    else
        looseinds=prot1+1:prot2-1;
    end
    % we need soft in between because they're on the same side of the 1st point
    for iPt=1:N
        indexm1=iPt-1;
        if (indexm1==0)
            indexm1=N;
        end 
        indexp1=iPt+1;
        if (indexp1==N+1)
            indexp1=1;
        end 
        if (iPt~=prot1 && iPt~=prot2)
            if (sum(looseinds==iPt) > 0) % iPt is in the soft section
                tao1=[x(iPt)-x(indexm1) y(iPt)-y(indexm1)]/...
                    norm([x0(iPt)-x0(indexm1) y0(iPt)-y0(indexm1)]);
                T1=gamsoft+ksoft*(norm(tao1)-1);
                tao1=tao1/norm(tao1);
                tao2=[x(indexp1)-x(iPt) y(indexp1)-y(iPt)]/...
                    norm([x0(iPt)-x0(indexp1) y0(iPt)-y0(indexp1)]);
                T2=gamsoft+ksoft*(norm(tao2)-1);
                tao2=tao2/norm(tao2);
                eforce(iPt,:)=(T2*tao2-T1*tao1)/h;
            else % iPt is in the hard section
                tao1=[x(iPt)-x(indexm1) y(iPt)-y(indexm1)]/...
                    norm([x0(iPt)-x0(indexm1) y0(iPt)-y0(indexm1)]);
                T1=gamhard+khard*(norm(tao1)-nrl);
                tao1=tao1/norm(tao1);
                tao2=[x(indexp1)-x(iPt) y(indexp1)-y(iPt)]/...
                    norm([x0(iPt)-x0(indexp1) y0(iPt)-y0(indexp1)]);
                T2=gamhard+khard*(norm(tao2)-nrl);
                tao2=tao2/norm(tao2);
                eforce(iPt,:)=(T2*tao2-T1*tao1)/h;
            end
        elseif (iPt==prot1) % it is prot1
            tao1=[x(iPt)-x(indexm1) y(iPt)-y(indexm1)]/...
                    norm([x0(iPt)-x0(indexm1) y0(iPt)-y0(indexm1)]);
            T1=gamhard+khard*(norm(tao1)-nrl); % hard behind it
            tao1=tao1/norm(tao1);
            tao2=[x(indexp1)-x(iPt) y(indexp1)-y(iPt)]/...
                    norm([x0(iPt)-x0(indexp1) y0(iPt)-y0(indexp1)]);
            T2=gamsoft+ksoft*(norm(tao2)-1); % soft ahead of it
            tao2=tao2/norm(tao2);
            eforce(iPt,:)=(T2*tao2-T1*tao1)/h;
        else % it is prot2
            tao1=[x(iPt)-x(indexm1) y(iPt)-y(indexm1)]/...
                    norm([x0(iPt)-x0(indexm1) y0(iPt)-y0(indexm1)]);
            T1=gamsoft+ksoft*(norm(tao1)-1); % soft behind it
            tao1=tao1/norm(tao1);
            tao2=[x(indexp1)-x(iPt) y(indexp1)-y(iPt)]/...
                    norm([x0(iPt)-x0(indexp1) y0(iPt)-y0(indexp1)]);
            T2=gamhard+khard*(norm(tao2)-nrl); % hard ahead of it
            tao2=tao2/norm(tao2);
            eforce(iPt,:)=(T2*tao2-T1*tao1)/h;
        end
    end
end 
        
    