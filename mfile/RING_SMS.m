function S = RING_SMS(kSpace, theta, Nsms)
% S = RING_SMS(kSpace,theta)
%    Inputs:
%       kSpace: radial k-space data
%          dimention: [sx, nor, coil]
%              sx: number of readouts along a ray
%              nor: total number of rays in temporal order
%              coil: number of coils (can be after coil compression)
%       theta: projection angle
%          dimention: [1, nor]
%
% Reference:
% Simple auto-calibrated gradient delay estimation from few spokes using
% Radial Intersections (RING)
% https://doi.org/10.1002/mrm.27506
%
% Modified for radial SMS data in CRIMP paper
% The code assumes Nsms = 3. Needs to be modified for more general use.
%
% Author: Ye Tian phye1988@gmail.com
%

N = 80; % number of rays used to calibration

for i=1:N
    % use the i ray, and find its 90 degree partener within 60 rays
    ray_pick = (i-1)*Nsms+1;
    ray_1 = squeeze(kSpace(:,ray_pick,:));
    theta_1 = theta(ray_pick);
    d_theta = theta_1 - theta(ray_pick-1+(Nsms+1:Nsms:60));
    [~,ray_pick_partner] = min(abs(abs(d_theta)-pi/2));
    ray_pick_partner = ray_pick_partner*Nsms+ray_pick;
    theta_1_partner = theta(ray_pick_partner);
    ray_1_partner = squeeze(kSpace(:,ray_pick_partner,:));
    
    % get the center part of picked k-space ray and its partener
    ray_1_im = fftshift(ifft(ray_1),1);
    ray_1_partner_im = fftshift(ifft(ray_1_partner),1);
    % remove outer part of image domain
    sx = size(ray_1_im,1);
    drop = round(sx*0.2);
    ray_1_im([1:drop,end-drop+1:end],:) = [];
    ray_1_partner_im([1:drop,end-drop+1:end],:) = [];
    % zero pad
    sx_new = size(ray_1_im,1);
    ray_1_im(end+1:sx*100,:) = 0;
    ray_1_im = circshift(ray_1_im,[round(sx*100/2),0]);
    ray_1_partner_im(end+1:sx*100,:) = 0;
    ray_1_partner_im = circshift(ray_1_partner_im,[round(sx*100/2),0]);
    
    % fft back to k-spac
    ray_1_k = fft(fftshift(ray_1_im,1));
    ray_1_partener_k = fft(fftshift(ray_1_partner_im,1));
    % keep center part
    center_new = round(sx*100/2)+1;
    ray_1_k = ray_1_k(round(center_new-sx_new*0.75):round(center_new+sx_new*0.75),:);
    ray_1_partener_k = ray_1_partener_k(round(center_new-sx_new*0.75):round(center_new+sx_new*0.75),:);
    % find the cross point
    d = ray_1_k-permute(ray_1_partener_k.',[3,1,2]);
    d = squeeze(sum(abs(d).^2,2));
    [~,x(i)] = min(min(d,[],2),[],1);
    [~,y(i)] = min(min(d,[],1),[],2);
    theta_pick(i) = theta_1;
    theta_partner(i) = theta_1_partner;
end

center = size(ray_1_k,1);
center = (center+1)/2;
x = (x-center)/100;%/(sx-round(sx*0.2)*2)*sx;
y = (y-center)/100;%/(sx-round(sx*0.2)*2)*sx;

e1 = cos(theta_pick) - cos(theta_partner);
e2 = sin(theta_pick) - sin(theta_partner);

A = zeros(2,3,N);
A(1,1,:) = e1;
A(1,3,:) = e2;
A(2,2,:) = e2;
A(2,3,:) = e1;

b = zeros(2,1,N);
b(1,1,:) = x.*cos(theta_pick) - y.*cos(theta_partner);
b(2,1,:) = x.*sin(theta_pick) - y.*sin(theta_partner);

A = permute(A,[1,3,2]);
A = reshape(A,2*N,3);
b = b(:);

s = A\b;

S = [s(1),s(3);s(3),s(2)];

end




