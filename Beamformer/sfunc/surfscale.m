function [Xout, Yout, Zout] = surfscale(q,X,Y,Z)
%surfscale - takes a surface scalar and surface vectors as input, returns
%input grid scaled by q as output.  Ensure X,Y,Z is 3D and has equal
%dimensions.

siz = size(X);
assert(length(siz) == 2,'Surface dimensions != 2');
assert(siz(1) == siz(2));

Xout = q.*X;
Yout = q.*Y;
Zout = q.*Z;

end