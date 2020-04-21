function Image = sr_forward(Image,iso)
[sx,sy,nof,Nslice] = size(Image);
Image = Image(:,:,:,iso.order_back);
Image = reshape(Image,[sx,sy,nof,1,iso.nSMS,iso.nset]);
end