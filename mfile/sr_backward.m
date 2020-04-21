function Image = sr_backward(Image,iso)

[sx,sy,nof,~,~,~] = size(Image);
Image = reshape(Image,[sx*sy*nof,iso.nSMS*iso.nset]);
Image = Image(:,iso.order);
Image = reshape(Image,[sx,sy,nof,iso.nSMS*iso.nset]);

end