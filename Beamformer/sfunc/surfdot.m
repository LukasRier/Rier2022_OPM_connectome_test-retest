function [D] = surfdot(q,X,Y,Z)
%surfdot - returns dot product of q with surface input at each point
%on the surface

siz = size(X);

assert(length(siz) == 2,'Surface dimensions != 2');
assert(siz(1) == siz(2));

D = zeros(siz);

%Iterate through all x, y, z values, make a vector v of each at each point,
%calculate q . v at each point, make a new meshgrid and populate with
%dprod values.
for ii = 1:siz(1)
    for jj = 1:siz(2)
            vect = [X(ii,jj); Y(ii,jj); Z(ii,jj)];
            dprod = dot(q,vect);
            D(ii,jj) = dprod;
    end
end