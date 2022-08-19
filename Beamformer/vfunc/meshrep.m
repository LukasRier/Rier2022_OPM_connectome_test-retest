function [Xout, Yout, Zout] = meshrep(q,X)
%meshrep - takes a single vector q and meshgrid vector as input, returns
%a grid of just q repeated at every point in space.  Ensure X,Y,Z is 3D and
%has equal dimensions.

siz = size(X);
assert(length(siz) == 3,'Meshgrid dimensions != 3');
assert(siz(1) == siz(2));
assert(siz(2) == siz(3));
assert(siz(3) == siz(1));

Xout = q(1)*ones(siz);
Yout = q(2)*ones(siz);
Zout = q(3)*ones(siz);

end