function [Bx_tot, By_tot, Bz_tot] = pointsBfield(Q,R0,Pts)

%Sanitization of input
sizQ = size(Q);
sizR0 = size(R0);
sizPts = size(Pts);

assert(sum(sizQ == sizR0) == length(sizQ), 'Dipole and position arrays do not match in size');
assert(length(sizQ) == 1 || length(sizQ) == 2);
assert(sizPts(2) == 3);

%Constants
mu_0 = 4*pi*1e-7;  %vacuum permeability
Qmag = 1e-9;  %|Q| = 5nAm

Bx_tot = zeros(sizPts(1),1);
By_tot = zeros(sizPts(1),1);
Bz_tot = zeros(sizPts(1),1);

xp = Pts(:,1);
yp = Pts(:,2);
zp = Pts(:,3);

R = sqrt(xp.^2 + yp.^2 + zp.^2);  %calculate |r| at each point

for n = 1:sizQ(1)

%Get points of a = r - R0
apx = xp - R0(n,1);
apy = yp - R0(n,2);
apz = zp - R0(n,3);

A = sqrt(apx.^2 + apy.^2 + apz.^2);  %calculate |a| = |r - Rq| at each point

R0x = repmat(R0(n,1),[sizPts(1) 1]);
R0y = repmat(R0(n,2),[sizPts(1) 1]);
R0z = repmat(R0(n,3),[sizPts(1) 1]);

% R0r = surfdot(R0(n,:),xsp,ysp,zsp);  %get r0 . r at each point

R0r = R0x.*xp + R0y.*yp + R0z.*zp;

% AR = surfdot2(xsp,ysp,zsp,aspx,aspy,aspz);

AR = xp.*apx + yp.*apy + zp.*apz;
F = A.*(R.*A + R.^2 - R0r);
rcoeff = (R.^-1).*(A.^2) + (A.^-1).*(AR) + 2*A + 2*R;
r0coeff = A + 2*R + (A.^-1).*AR;

rpx = xp.*rcoeff;
rpy = yp.*rcoeff;
rpz = zp.*rcoeff;

r0px = R0x.*r0coeff;
r0py = R0y.*r0coeff;
r0pz = R0z.*r0coeff;

dFx = rpx - r0px;  dFy = rpy - r0py;  dFz = rpz - r0pz;

cqr0 = repmat(cross(Q(n,:),R0(n,:)), [sizPts(1) 1]);
cqr0x = cqr0(:,1);
cqr0y = cqr0(:,2);
cqr0z = cqr0(:,3);

% [cqr0x, cqr0y, cqr0z] = surfrep(cross(Q(n,:),R0(n,:)),xsp);
% [Fcqr0x, Fcqr0y, Fcqr0z] = surfscale(F,cqr0x,cqr0y,cqr0z);

Fcqr0x = F.*cqr0x;
Fcqr0y = F.*cqr0y;
Fcqr0z = F.*cqr0z;

% qr0r = surfdot(cross(Q(n,:),R0(n,:)),xsp,ysp,zsp);

qr0r = cqr0x.*xp + cqr0y.*yp + cqr0z.*zp;

%Full B-field calculation
Bx = (F.^-2).*(Fcqr0x - qr0r.*dFx);
By = (F.^-2).*(Fcqr0y - qr0r.*dFy);
Bz = (F.^-2).*(Fcqr0z - qr0r.*dFz);

% size(Bx)
% size(Bx_tot)

Bx_tot = Bx_tot + Bx;
By_tot = By_tot + By;
Bz_tot = Bz_tot + Bz;

end

Bx_tot = (Qmag*mu_0/(4*pi))*Bx_tot;
By_tot = (Qmag*mu_0/(4*pi))*By_tot;
Bz_tot = (Qmag*mu_0/(4*pi))*Bz_tot;



end