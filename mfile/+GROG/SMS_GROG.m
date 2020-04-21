function [kSpace_cart,G] = SMS_GROG(kSpace_radial, kx, ky, phase_mod, para, G_pre_calculated)
%  This is a GROG interpolation function for 2D MRI k-space data
%
%  Reference:
%  Ye Tian, et al. (2019) PLOS ONE 14(2): e0211738.
%  Nicole Seiberlich, et al. (2008) MRM 59:930-935.
%  
%  Inputs:
%         kSpace_radial: up to 5 dimentions
%         kx, ky       : k-space positions
%
% This version does not support multi-slices inputs.
%
% Copyright: University of Utah Cardiovascular MRI Group
% https://medicine.utah.edu/radiology/radiology-research/research-labs/dibella2/
%
% Contact:
% Ye Tian, phye1988@gmail.com




nSMS = max(phase_mod(:))+1;

[sx,nor,nof,nc,~,NSlice] = size(kSpace_radial);
core_size = para.core_size;
over_sampling = para.over_sampling;

G = cell(1,nSMS);

for j=1:nSMS
    for i=1:NSlice
        kSpace_radial_temp = kSpace_radial(:,:,:,:,i);
        SMS_subset = phase_mod == j-1;
        kSpace_radial_temp = kSpace_radial_temp(repmat(SMS_subset,[sx 1 1 nc]));
        kSpace_radial_temp = reshape(kSpace_radial_temp,[sx,nor/nSMS,nof,nc]);
        kx_temp = kx(repmat(SMS_subset,[sx 1 1]));
        ky_temp = ky(repmat(SMS_subset,[sx 1 1]));
        kx_temp = reshape(kx_temp,[sx nor/nSMS nof]);
        ky_temp = reshape(ky_temp,[sx nor/nSMS nof]);

        if exist('G_pre_calculated')
            G{j} = G_pre_calculated{j};
            kSpace_cart(:,:,:,:,:,:,j) = GROG.GROG_Dictionary_interp_new(kSpace_radial_temp,G_pre_calculated{j}.Gx,G_pre_calculated{j}.Gy,kx_temp,ky_temp);
        else
            G{j} = GROG.GNUFFT_init(kSpace_radial_temp,kx_temp,ky_temp,over_sampling,core_size);
            kSpace_cart(:,:,:,:,:,:,j) = GROG.GNUFFT_rad2cart(kSpace_radial_temp,G{j});
           
        end

    end
end

sx_over = size(kSpace_cart,1);
kSpace_cart = reshape(kSpace_cart,[sx_over,sx_over,nof,nc,1,NSlice,nSMS]);
if size(kSpace_cart,1)~=sx*over_sampling
    kSpace_cart([1,end],:,:,:,:,:,:) = [];
    kSpace_cart(:,[1,end],:,:,:,:,:) = [];
end
