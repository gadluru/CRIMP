function Image = sr_forward(Image,siz)
[sx,sy,nof,~] = size(Image);
Image = Image(:,:,:,siz.order_back);
Image = reshape(Image,[sx,sy,nof,1,siz.nSMS,siz.nset]);
end