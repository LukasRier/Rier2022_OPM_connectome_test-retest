function [D] = meshdot(q,X,Y,Z)
%meshdot - returns dot product of q with meshgrid input at each point
%on the grid.

siz = size(X);

assert(length(siz) == 3,'Meshgrid dimensions != 3');
assert(siz(1) == siz(2));
assert(siz(2) == siz(3));
assert(siz(3) == siz(1));

D = zeros(siz);

%Iterate through all x, y, z values, make a vector v of each at each point,
%calculate q . v at each point, make a new meshgrid and populate with
%dprod values.
for zn = 1:siz(1)
    for yn = 1:siz(1)
        for xn = 1:siz(1)
            vect = [X(xn,yn,zn); Y(xn,yn,zn); Z(xn,yn,zn)];
            dprod = dot(q,vect);
            D(xn,yn,zn) = dprod;
            
        end
    end
end