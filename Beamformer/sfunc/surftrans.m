function [Xout, Yout, Zout] = surftrans(q,X,Y,Z)
%meshtrans - takes a single vector q and surface arrays as input, returns
%input surface translated by q as output.  Ensure X,Y,Z are 2D and have equal
%dimensions.

siz = size(X);
assert(length(siz) == 2,'Surface dimensions != 2');
assert(siz(1) == siz(2));

Xout = X + q(1)*ones(siz);
Yout = Y + q(2)*ones(siz);
Zout = Z + q(3)*ones(siz);

end