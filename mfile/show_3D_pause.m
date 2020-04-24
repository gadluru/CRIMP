function show_3D_pause(img, brightness)
% show_3D_yt_pause(img)
% display 3D dynamic images frame-by-frame

if nargin==1
    brightness = 0.2;
end
[sx,sy,ns,nf] = size(img);
ns_sqrt = ceil(sqrt(ns));
if ns_sqrt^2~=ns
    img(:,:,ns+1:ns_sqrt^2,:) = zeros(sx,sy,ns_sqrt^2-ns,nf);
    ns = ns_sqrt^2;
end
img = reshape(img,sx,sy*ns_sqrt,ns_sqrt,nf);
img = permute(img,[1 3 2 4]);
img = reshape(img,sx*ns_sqrt,sy*ns_sqrt,nf);
figure
for i=1:size(img,3)
    imagesc(abs(img(:,:,i)))
    colormap gray
    axis image
    brighten(brightness)
    title(i)
    %drawnow
    pause
end
end