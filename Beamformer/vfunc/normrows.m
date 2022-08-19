function [out] = normrows(in)

assert(length(size(in)) == 2);

sq2 = in.^2;
nr2 = sqrt(sum(sq2,2));
out = in./repmat(nr2,[1 size(in,2)]);

end