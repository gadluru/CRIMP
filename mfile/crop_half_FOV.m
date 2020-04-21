function Image = crop_half_FOV(Image)

sx = size(Image,1);
sy = size(Image,2);
cut_x = round(sx/4);
cut_y = round(sy/4);
Image(1:cut_x,:,:,:,:,:,:) = [];
Image(cut_x*2+1:end,:,:,:,:,:,:) = [];
Image(:,1:cut_y,:,:,:,:,:) = [];
Image(:,cut_y*2+1:end,:,:,:,:,:) = [];
