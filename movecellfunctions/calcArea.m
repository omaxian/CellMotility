% Use Green's theorem to calculate the area of the cell from the x
% and y locations.
function area = calcArea(x,y)
    N=length(x);
    area=0;
    for iPt=1:N-1
        area=area+0.5*(x(iPt)*y(iPt+1)-y(iPt)*x(iPt+1));
    end
    area=area+0.5*(x(N)*y(1)-y(N)*x(1));
end 