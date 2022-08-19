function [ outmatrix ] = rotateAngleMatrix( u, angle )
%ROTATEANGLEMATRIX Summary of this function goes here
%   Detailed explanation goes here

cprodmat = [ 0 -u(3) u(2); ...
            u(3) 0 -u(1);  ...
            -u(2) u(1) 0];

kronu = kron(u,u');

R = cosd(angle) * eye(3) + sind(angle) * cprodmat + (1 - cosd(angle))* kronu;

outmatrix = eye(4);
outmatrix(1:3, 1:3) = R;



end

