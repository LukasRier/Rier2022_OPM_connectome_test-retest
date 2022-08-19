function [vtx, vty, vtz] = dipolefangrid(rx,ry,rz,theta_t)

sizr = size(rx);

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

sizt = [sizr ntheta];
vtx = zeros(sizt);
vty = zeros(sizt);
vtz = zeros(sizt);

for zi = 1:sizr(1)
    for yi = 1:sizr(2)
        for xi = 1:sizr(3)
% thx1 = thx(xi,yi,zi);
% thy1 = thy(xi,yi,zi);
% thz1 = thz(xi,yi,zi);
% phx1 = phx(xi,yi,zi);
% phy1 = phy(xi,yi,zi);
% phz1 = phz(xi,yi,zi);

thx1 = thx(xi,yi,zi);
thy1 = thy(xi,yi,zi);
thz1 = thz(xi,yi,zi);
phx1 = phx(xi,yi,zi);
phy1 = phy(xi,yi,zi);
phz1 = phz(xi,yi,zi);


vtx(xi,yi,zi,:) = cos(theta_t).*thx1 + sin(theta_t).*phx1;
vty(xi,yi,zi,:) = cos(theta_t).*thy1 + sin(theta_t).*phy1;
vtz(xi,yi,zi,:) = cos(theta_t).*thz1 + sin(theta_t).*phz1;

% vtx(xi,yi,zi,:) = thx1;
% vty(xi,yi,zi,:) = thy1;
% vtz(xi,yi,zi,:) = thz1;
        
        end
    end
end

end