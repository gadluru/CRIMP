function showImage3D(Image,Cost)
Image = crop_half_FOV(abs(Image));
[nx,ny,nof,nz] = size(Image);
frame_num = floor(nof/4);

im = Image(:,:,[frame_num, frame_num*2, frame_num*3],round(nz/2));
im = permute(im,[1 3 2]);
im = reshape(im,[nx*3 ny]);

figure(100)
subplot(1,2,1)
imagesc(im)
colormap gray
brighten(0.3)
axis image
axis off
        
subplot(1,2,2)
plotCost(Cost)
drawnow
end
