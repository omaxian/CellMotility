% Routine to fix any geometrical issues that result from a coarse
% discretization, in particular if a cortex node is inside the nucleus or
% an ECM node inside the cortex. 
% Inputs: the points that we are checking, (xpoly,ypoly) = vertices of the
% shape we are checking if points are inside, np1in = initial guess for the
% normal vector, cordist = distance we make random moves. 
function newlocs = correctPoints(points,xpoly,ypoly,np1in,cordist)
    np1=np1in;
    [in,edg]=inpolygon(points(:,1),points(:,2),[xpoly; xpoly(1)],[ypoly; ypoly(1)]);
    while (sum(in-edg) > 0) % if it's inside, keep picking random vector until outside
         inNode = find((in-edg)==1);
         points(inNode,:)=points(inNode,:)-cordist*np1;
         np1 = rand(1,2)-0.5;
         np1 = np1/norm(np1);
         points(inNode,:)=points(inNode,:)+cordist*np1;
         [in,edg]=inpolygon(points(:,1),points(:,2),[xpoly; xpoly(1)],[ypoly; ypoly(1)]);
    end
    newlocs = points;
end
    