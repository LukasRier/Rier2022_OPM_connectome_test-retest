function [xtc, ytc, ztc] = findxyz(rx, ry, rz, dlocx, dlocy, dlocz)

sizr = size(rx);

[~,Iy] = min(abs(rx(:) - dlocx));
[~,Ix] = min(abs(ry(:) - dlocy));
[~,Iz] = min(abs(rz(:) - dlocz));

[xi, xj, xk] = ind2sub(sizr,Ix);
[yi, yj, yk] = ind2sub(sizr,Iy);
[zi, zj, zk] = ind2sub(sizr,Iz);

xtc = max([xi, xj, xk]);
ytc = max([yi, yj, yk]);
ztc = max([zi, zj, zk]);

end