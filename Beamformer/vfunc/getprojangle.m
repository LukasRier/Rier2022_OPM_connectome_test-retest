function [ theta_t ] = getprojangle( Q, R0)
%find angle of projection of dipole onto theta-phi plane

theta_t = zeros(size(Q,1));

R = sqrt(R0(:,1).^2 + R0(:,2).^2 + R0(:,3).^2);
theta = acos(R0(:,3)./R);
phi = atan2(R0(:,2),R0(:,1));

%Get theta unit vector
thx = cos(theta).*cos(phi);
thy = cos(theta).*sin(phi);
thz = -sin(theta);

%Get phi unit vector
phx = -sin(phi);
phy = cos(phi);
phz = zeros(size(phx));

u = dot(Q,[thx, thy, thz],2);
v = dot(Q,[phx, phy, phz],2);

for n = 1:size(u,1)
theta_t(n) = atan2d(u(n),v(n));
end

end

