function [Image_reg,shifts,mask] = rigid_reg_3D(Image,mask)

[sx,sy,sz,nof] = size(Image);
if ~exist('mask')
    for i=1:sz
        mask(:,:,i) = get_mask(Image(:,:,i,end));
    end
end


%Image = cat(3,zeros(sx,sy,1,nof),Image,zeros(sx,sy,1,nof));
Image = cat(3,Image(:,:,[1,1],:),Image,Image(:,:,[end,end,end],:));
Image = interp_2_slices_in_between(Image);
Image_reg = Image;

mask = cat(3,zeros(sx,sy,2),mask,zeros(sx,sy,3));
mask = interp_2_slices_in_between(mask);
mask(:,:,[1:2,end-2:end]) = false;

sz = size(Image,3);
noi = 4;
step = 2;
shifts = zeros(3,noi,nof);
first_iter_flag = 1;

step = ones(1,nof)*step;
z = round(sz/2+1);
[x,y] = find(mask(:,:,z));
x = round(mean(x));
y = round(mean(y));

for iter_no = 1:noi
    Image_ref = Image_reg(:,:,:,end);
    for i=nof-1:-1:1
        tic;
        [Image_reg(:,:,:,i),shifts(:,iter_no,i),step(i)] = Reg_rigid_3D_fast(Image(:,:,:,i),Image_ref,step(i),10,mask,first_iter_flag);
        Image_ref = (Image_ref*(nof-i) + Image_reg(:,:,:,i))/(nof-i+1);
        %mean(Image_reg(:,:,:,i+1:end),4)
%                 figure(1)
%                 line_x = squeeze(Image_reg(:,y,z,:));
%                 line_y = squeeze(Image_reg(x,:,z,:));
%                 line_z = squeeze(Image_reg(x,y,:,:));
%                 subplot(5,1,3)
%                 imagesc(line_x);
%                 subplot(5,1,4)
%                 imagesc(line_y);
%                 subplot(5,1,5)
%                 imagesc(line_z);
%                 drawnow
        fprintf(['Iteration = ' num2str(iter_no) ', ' 'Frame = ' num2str(i) '...']);
        toc;
%         figure(1)
%         imagesc(mean(Image_reg(:,:,z,i+1:end),4))
%         axis image
%         colormap gray
%         brighten(0.4)
%         drawnow


    end
    Image = Image_reg;
%     step = step/2;
    first_iter_flag = 0;

end

Image_reg = Image_reg(:,:,1:3:end,:);
mask = mask(:,:,1:3:end);
Image_reg = Image_reg(:,:,3:end-3,:);
mask = mask(:,:,3:end-3,:);

end


function Image = interp_2_slices_in_between(Image)
% Image should have dimention: x,y,z,frames
[sx,sy,sz,nof] = size(Image);
if nof>1
    [x,y,z,f] = ndgrid(1:sx,1:sy,1:sz,1:nof);
    [xv,yv,zv,fv] = ndgrid(1:sx,1:sy,1:1/3:sz,1:nof);
    Image = interpn(x,y,z,f,Image,xv,yv,zv,fv);
else
    [x,y,z] = ndgrid(1:sx,1:sy,1:sz);
    [xv,yv,zv] = ndgrid(1:sx,1:sy,1:1/3:sz);
    Image = interpn(x,y,z,Image,xv,yv,zv);
end
end