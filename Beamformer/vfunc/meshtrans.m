function [Xout, Yout, Zout] = meshtrans(q,X,Y,Z)
%meshtrans - takes a single vector q and meshgrid vectors as input, returns
%input grid translated by q as output.  Ensure X,Y,Z is 3D and has equal
%dimensions.

siz = size(X);
assert(length(siz) == 3,'Meshgrid dimensions != 3');
assert(siz(1) == siz(2));
assert(siz(2) == siz(3));
assert(siz(3) == siz(1));

Xout = X + q(1)*ones(siz);
Yout = Y + q(2)*ones(siz);
Zout = Z + q(3)*ones(siz);

end