function [vtx, vty, vtz] = tangdipoles(rx,ry,rz,theta_t)

sizr = size(rx);

sizt = size(theta_t);
assert(sum(sizt) == 2);

% GRIDS
% Get back R and unit vectors from passed grid
R = sqrt(rx.^2 + ry.^2 + rz.^2);  %calculate |r| at each point

%getting normal vectors
erx = rx./R;
ery = ry./R;
erz = rz./R;

%Get theta and phi for each point in grid
theta = acos(rz./R);
phi = atan2(ry,rx);

%Get theta unit vector at every point
thx = cos(theta).*cos(phi);
thy = cos(theta).*sin(phi);
thz = -sin(theta);

%Get phi unit vector at every point
phx = -sin(phi);
phy = cos(phi);
phz = zeros(size(rx));

ntheta = size(theta_t,2);
% theta_t = linspace(0, pi, ntheta);

vtx = zeros(sizr);
vty = zeros(sizr);
vtz = zeros(sizr);


vtx = cos(theta_t).*thx + sin(theta_t).*phx;
vty = cos(theta_t).*thy + sin(theta_t).*phy;
vtz = cos(theta_t).*thz + sin(theta_t).*phz;

end