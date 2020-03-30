% Add the pinch to the pulling force.
% Uses reference normal vector at jPt for the protrusive force
% Inputs: fmag = f0 in the paper, magnitude of force, N = number of points,
% s = arclength parameterization, jPt = the point on the cortex where
% protrusion forms
% Outputs: the new force
function extraforce=addpinch(fmag,N,s,jPt)
    extraforce=zeros(N,2);
    normals=[cos(s) sin(s)];
    indices=[jPt-1:jPt+1];
    indices=mod(indices,N);
    indices(indices==0)=N;
    normals=normals(indices,:)';
    otherpts=setdiff((1:N),indices);
    values=[0.5 1 0.5]';
    extraforce(indices,:)=fmag*values.*normals';
    extraforce(otherpts,:)=-sum(extraforce)/(N-3)+ extraforce(otherpts,:);
end