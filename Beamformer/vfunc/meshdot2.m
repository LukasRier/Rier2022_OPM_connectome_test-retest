function [D] = meshdot2(X,Y,Z,U,V,W)
%meshdot2 - returns dot product grid of two vector grids.

siz = size(X);
siz2 = size(U);

assert(length(siz) == 3,'Meshgrid 1 dimensions != 3');
assert(length(siz2) == 3,'Meshgrid 2 dimensions != 3');
assert(siz(1) == siz(2));
assert(siz(2) == siz(3));
assert(siz(3) == siz(1));

D = zeros(siz);

dotX = X.*U;
dotY = Y.*V;
dotZ = Z.*W;

D = dotX + dotY + dotZ;