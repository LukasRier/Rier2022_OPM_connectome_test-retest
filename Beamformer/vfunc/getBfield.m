function [Bx_tot, By_tot, Bz_tot] = getBfield(Q,R0,rx,ry,rz)

%Sanitization of input
sizQ = size(Q);
sizR0 = size(R0);
sizrx = size(rx);
sizry = size(ry);
% sizrz = size(rz);

assert(sum(sizQ == sizR0) == length(sizQ), 'Dipole and position arrays do not match in size');
assert(sum(sizrx == sizry) == length(sizrx), 'Supplied grid arrays do not match in size');
assert(length(sizQ) == 1 || length(sizQ) == 2);

%Constants
mu_0 = 4*pi*1e-7;  %vacuum permeability
Qmag = 1e-9;  %|Q| = 1nAm



%GRIDS
%Create r and start grid centred on 0
R = sqrt(rx.^2 + ry.^2 + rz.^2);  %calculate |r| at each point
Bx_tot = zeros(size(rx));
By_tot = zeros(size(rx));
Bz_tot = zeros(size(rx));

for n = 1:sizQ(1)
    
%Get grid of a = r - R0
[ax,ay,az] = meshtrans(-R0(n,:),rx,ry,rz);
A = sqrt(ax.^2 + ay.^2 + az.^2);  %calculate |a| = |r - Rq| at each point
[R0x, R0y, R0z] = meshrep(R0(n,:),rx); %make R0 field
R0r = meshdot(R0(n,:),rx,ry,rz);  %get r0 . r at each point
AR = meshdot2(rx,ry,rz,ax,ay,az);
F = A.*(R.*A + R.^2 - R0r);
rcoeff = (R.^-1).*(A.^2) + (A.^-1).*(AR) + 2*A + 2*R;
r0coeff = A + 2*R + (A.^-1).*AR;
[rpx, rpy, rpz] = meshscale(rcoeff,rx,ry,rz);
[r0px, r0py, r0pz] = meshscale(r0coeff,R0x,R0y,R0z);
dFx = rpx - r0px;  dFy = rpy - r0py;  dFz = rpz - r0pz;
[cqr0x, cqr0y, cqr0z] = meshrep(cross(Q(n,:),R0(n,:)),rx);
[Fcqr0x, Fcqr0y, Fcqr0z] = meshscale(F,cqr0x,cqr0y,cqr0z);

qr0r = meshdot(cross(Q(n,:),R0(n,:)),rx,ry,rz);

%Full B-field calculation
Bx = (F.^-2).*(Fcqr0x - qr0r.*dFx);
By = (F.^-2).*(Fcqr0y - qr0r.*dFy);
Bz = (F.^-2).*(Fcqr0z - qr0r.*dFz);

Bx_tot = Bx_tot + Bx;
By_tot = By_tot + By;
Bz_tot = Bz_tot + Bz;

    
end

Bx_tot = (Qmag*mu_0/(4*pi))*Bx_tot;
By_tot = (Qmag*mu_0/(4*pi))*By_tot;
Bz_tot = (Qmag*mu_0/(4*pi))*Bz_tot;



end