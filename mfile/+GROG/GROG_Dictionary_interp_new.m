function kSpace_cart = GROG_Dictionary_interp_new(kSpace_radial,Gx,Gy,kx,ky)
% kSpace_cart = GROG_Dictionary_interp_new(kSpace_radial,Gx,Gy,kx,ky)
% Pre-calculate a dictionary to store all the possible shifts for GROG
% operator from -0.50 to 0.50 at 0.01 intervel. A lot faster than calculate
% the shift operator every time.
%
% inputs:
%   kSpace_radial: input k space data, should have size 
%       dimentions: [sx, nor, nof, nc]
%   Gx, Gy: pre-calculated GROG operator in x and y directions
%   kx, ky: k space trajectory, have value [-sx/2:sx/2]

% remind that this can also handle asymetric k space data, however the kx
% ky shoud be full length. This code uses size of kx ky to be the size of
% cart space.

% copyright Ye Tian phye1988@gmail.com
% University of Utah DiBella Group

[sx,nor,nof,nc] = size(kSpace_radial);
kSpace_radial = reshape(kSpace_radial,[sx*nor*nof,1,nc]);
% in case of asym kspace. skx shoule be the cart size that interpolate on
skx = size(kx,1);

% pre-allocate kSpace and weight in Cart space
%kSpace_cart = zeros([(skx+3)*(skx+3) nof nc],'single');
%weight_cart = zeros([(skx+3)*(skx+3) nof],'single');
%if isa(kSpace_radial,'gpuArray')
%    kSpace_cart = gpuArray(kSpace_cart);
%    weight_cart = gpuArray(weight_cart);
%end

% get shift distance
kx_round = round(kx);
ky_round = round(ky);
dx = round((kx_round - kx)*100)/100;
dy = round((ky_round - ky)*100)/100;
weight = 1 - sqrt(dx.^2 + dy.^2);

% calculate the dictionary for GROG operator
GxDict = single(zeros([nc nc 101]));
GyDict = single(zeros([nc nc 101]));
if isa(Gx,'gpuArray')
    GxDict = gpuArray(GxDict);
    GyDict = gpuArray(GyDict);
end

d = -0.50:0.01:0.50;
for di = 1:101
    GxDict(:,:,di) = Gx^d(di);
    GyDict(:,:,di) = Gy^d(di);
end

% calculate shift operator for all radial points
x_cart = kx_round + skx/2 + 2;
y_cart = ky_round + skx/2 + 2;

xDict = round((dx + 0.5)*100+1);
yDict = round((dy + 0.5)*100+1);

xDict = xDict(:);
yDict = yDict(:);

Gx_shift_all = GxDict(:,:,xDict);
Gx_shift_all = permute(Gx_shift_all,[3 1 2]);
Gy_shift_all = GyDict(:,:,yDict);
Gy_shift_all = permute(Gy_shift_all,[3 4 1 2]);
G_shift = bsxfun(@times,Gx_shift_all,Gy_shift_all);
G_shift = squeeze(sum(G_shift,3));

% calculate shifted radial points
k_target = bsxfun(@times,G_shift,kSpace_radial); clear kSpace_radial G*
k_target = squeeze(sum(k_target,3));
k_target = bsxfun(@times,weight(:),k_target);
k_target = reshape(k_target,[sx nor nof nc]);

indx = sub2ind([(skx+3),(skx+3),nof,1],x_cart,y_cart);
%k_target = permute(k_target,[1 2 4 3]);

%%%
k_target = reshape(k_target,[sx*nor*nof,nc]);
weight = reshape(weight,skx*nor,nof);
indx = reshape(indx,[skx*nor,nof]);

indx = bsxfun(@plus,indx,0:(skx+3)^2:(skx+3)^2*(nof-1));
radial_num = bsxfun(@plus,(1:skx*nor).',0:skx*nor:skx*nor*(nof-1));
rad2cart = sparse(indx(:),radial_num(:),1,(skx+3)^2*nof,skx*nor*nof);
kSpace_cart = single(rad2cart*double(k_target));
weight_cart = single(rad2cart*double(weight(:)));
%rad2cart = cell(1,nof);

%for i=1:nof
%    rad2cart_temp = sparse(indx(:,i),(1:skx*nor).',1,(skx+3)^2,skx*nor);
%    kSpace_cart(:,i,:) = single(rad2cart_temp*double(k_target(:,:,i)));
%    weight_cart(:,i) = single(rad2cart_temp*double(weight(:,i)));
%end

% apply the weight and cut the end of cart space
weight_cart(weight_cart==0)=1;
kSpace_cart = bsxfun(@rdivide,kSpace_cart,weight_cart);
kSpace_cart = reshape(kSpace_cart,[(skx+3) (skx+3) nof nc]);
kSpace_cart(1,:,:,:,:,:,:) = [];
kSpace_cart(:,1,:,:,:,:,:) = [];
kSpace_cart(end-1:end,:,:,:,:,:,:) = [];
kSpace_cart(:,end-1:end,:,:,:,:,:) = [];