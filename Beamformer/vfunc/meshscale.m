function [Xout, Yout, Zout] = meshscale(q,X,Y,Z)
%meshscale - takes a scalar grid q and meshgrid vectors as input, returns
%input grid scaled by q as output.  Ensure X,Y,Z is 3D and has equal
%dimensions.

siz = size(X);
assert(length(siz) == 3,'Meshgrid dimensions != 3');
assert(siz(1) == siz(2));
assert(siz(2) == siz(3));
assert(siz(3) == siz(1));

Xout = q.*X;
Yout = q.*Y;
Zout = q.*Z;

end