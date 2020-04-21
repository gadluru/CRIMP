function Image_lr = LowRank_patch_yt(Image)

siz = size(Image);
% if length(siz) < 4
%     siz(4) = 1;
% end

patch_size = [5,5,2];
patch_begin_x = 1:patch_size(1):siz(1);
if patch_begin_x(end) + patch_size(1) -1 > siz(1)
    patch_begin_x(end) = siz(1) - patch_size(1) + 1;
end
patch_begin_y = 1:patch_size(2):siz(2);
if patch_begin_y(end) + patch_size(2) -1 > siz(2)
    patch_begin_y(end) = siz(2) - patch_size(2) + 1;
end
patch_begin_z = 1:patch_size(3):siz(3);
if patch_begin_z(end) + patch_size(3) -1 > siz(3)
    patch_begin_z(end) = siz(3) - patch_size(3) + 1;
end

Image_lr = zeros(size(Image));

for nphase = 1:siz(5)
    Image_phase = Image(:,:,:,:,nphase);
    Image_lr_phase = zeros(siz(1:4),'single');
    Mask = Image_lr_phase;
    
    for nx = patch_begin_x
        for ny = patch_begin_y
            for nz = patch_begin_z
                patch_temp = Image_phase(nx:nx+patch_size(1)-1,ny:ny+patch_size(2)-1,nz:nz+patch_size(3)-1,:);
                patch_temp = reshape(patch_temp,[prod(patch_size),siz(4)]);
                [U,S,V] = svd(patch_temp,0);
                S = S - S(1)*0.02;
                S(S<0) = 0;
                patch_temp = U*S*V';
                patch_temp = reshape(patch_temp,[patch_size,siz(4)]);
                Image_lr_phase(nx:nx+patch_size(1)-1,ny:ny+patch_size(2)-1,nz:nz+patch_size(3)-1,:) = Image_lr_phase(nx:nx+patch_size(1)-1,ny:ny+patch_size(2)-1,nz:nz+patch_size(3)-1,:) + patch_temp;
                Mask(nx:nx+patch_size(1)-1,ny:ny+patch_size(2)-1,nz:nz+patch_size(3)-1,:) = Mask(nx:nx+patch_size(1)-1,ny:ny+patch_size(2)-1,nz:nz+patch_size(3)-1,:) + 1;
            end
        end
    end
    Mask = 1./Mask;
    Image_lr_phase = Image_lr_phase.*Mask;
    Image_lr(:,:,:,:,nphase) = Image_lr_phase; 
end

end
