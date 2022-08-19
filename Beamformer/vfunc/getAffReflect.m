function [ tform] = getAffReflect( plane_normal )
%GETAFFREFLECT Summary of this function goes here
%   Detailed explanation goes here

plane_normal = plane_normal(:);

reflmat = eye(4);
reflmat(1:3,1:3) = eye(3) - (2 * plane_normal * plane_normal');


tform = affine3d(reflmat);

end

