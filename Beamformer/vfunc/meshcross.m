function [Xcr, Ycr, Zcr] = meshcross(q,X,Y,Z)
%meshcross - returns cross product of q with meshgrid input at each point
%on the grid.

siz = size(X);

assert(length(siz) == 3,'Meshgrid dimensions != 3');
assert(siz(1) == siz(2));
assert(siz(2) == siz(3));
assert(siz(3) == siz(1));

%Iterate through all x, y, z values, make a vector v of each at each point,
%calculate q x v at each point, make a new meshgrid and populate each with
%x value, y values, z values.
for zn = 1:siz(1)
    for yn = 1:siz(1)
        for xn = 1:siz(1)
            vect = [X(xn,yn,zn); Y(xn,yn,zn); Z(xn,yn,zn)];
            cprod = cross(q,vect);
            Xcr(xn,yn,zn) = cprod(1);
            Ycr(xn,yn,zn) = cprod(2);
            Zcr(xn,yn,zn) = cprod(3);
            
        end
    end
end

end