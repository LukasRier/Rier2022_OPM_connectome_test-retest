function [D] = meshdot2(X,Y,Z,U,V,W)
%surfdot2 - returns dot product surface of two surfaces.

siz = size(X);
siz2 = size(U);

assert(length(siz) == 2,'Surface 1 dimensions != 2');
assert(length(siz2) == 2,'Surface 2 dimensions != 2');
assert(siz(1) == siz(2));
assert(siz2(2) == siz2(1));
assert(siz(1) == siz2(1));

D = zeros(siz);

dotX = X.*U;
dotY = Y.*V;
dotZ = Z.*W;

D = dotX + dotY + dotZ;