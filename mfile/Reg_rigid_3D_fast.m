function [I2,shifts,step] = Reg_rigid_3D_fast(I1,I0,step,noi,mask,flag)
[sx,sy,sz] = size(I0);

a = 30;
s = 3;

% for i=1:sz
%     mask(:,:,i) = get_mask(I0(:,:,i));
% end


kx = cos(2*pi*(0:sx-1)/sx);
ky = cos(2*pi*(0:sy-1)/sy);
kz = cos(2*pi*(0:sz/(sz-1):sz));
W = 2*(kx+ky.'+permute(kz,[3 1 2])-3);
W = (1-a*W).^-s;

N = sum(mask(:));
%[y,x,z] = meshgrid(1:sx,1:sy,1:sz);

I2_fft = fftshift3(fft3(fftshift3(I1)));

PhaseX = reshape(-2i*pi*((1:sx)-floor(sx/2)-1)/sx,sx,1);
PhaseY = reshape(-2i*pi*((1:sy)-floor(sy/2)-1)/sy,1,sy);
PhaseZ = reshape(-2i*pi*((1:sz)-floor(sz/2)-1)/sz,1,1,sz);

Phase = exp(PhaseX*0 + PhaseY*0 + PhaseZ*0);
shifts = zeros(3,1);
for iter=1:noi

%     I2 = interp3(I1,y,x,z);
    I2 = ifftshift3(ifft3(fftshift3(I2_fft.*Phase)));
    I2 = abs(I2);

%     figure(1)
%     imagesc(I2(:,:,5))
%     colormap gray
%     brighten(0.4)
%     axis image
%     drawnow
    
    ddx = 0.5*(I2(3:end,:,:) - I2(1:end-2,:,:));
    ddy = 0.5*(I2(:,3:end,:) - I2(:,1:end-2,:));
    ddz = 0.5*(I2(:,:,3:end) - I2(:,:,1:end-2));
    ddx = cat(1,I2(2,:,:) - I2(1,:,:),ddx,I2(end,:,:) - I2(end-1,:,:));
    ddy = cat(2,I2(:,2,:) - I2(:,1,:),ddy,I2(:,end,:) - I2(:,end-1,:));
    ddz = cat(3,I2(:,:,2) - I2(:,:,1),ddz,I2(:,:,end) - I2(:,:,end-1));
    
    dI = I2-I0;
    %dI = dI - medfilt3(dI,[3,3,3]);
    
    dx = W.*fft3(ddx.*dI);
    dx = -real(ifft3(dx));

    dy = W.*fft3(ddy.*dI);
    dy = -real(ifft3(dy));
    
    dz = W.*fft3(ddz.*dI);
    dz = -real(ifft3(dz));
    
    dx = dx.*mask;
    dy = dy.*mask;
    dz = dz.*mask;
    
%     d = sqrt(dx.^2+dy.^2+dz.^2);
%     md = max(d(:));
%     
%     dx = dx/md;
%     dy = dy/md;
%     dz = dz/md;


    
    
%     figure(2)
%     plot(E)
%     drawnow

    dx = -mean(dx(:));
    dy = -mean(dy(:));
    dz = -mean(dz(:))/4;

    d = sos([dx,dy,dz]);
    if d<0.01
        return
    end
    if iter==1 && flag
        step = step/d;
    end
    E(iter) = sum(d(:));
    if iter>1 && E(iter)>E(iter-1)
        step = step/2;
        %return
    end
    
    dx = dx*step;
    dy = dy*step;
    dz = dz*step;
    
    
    
    shifts = shifts + [dx;dy;dz];
    Phase = exp(PhaseX*shifts(1) + PhaseY*shifts(2) + PhaseZ*shifts(3));

%     I2_fft = fftshift3(fft3(fftshift3(I2)));
%     x = x + step*dx;
%     y = y + step*dy;
%     z = z + step*dz;
%     
%     x(x<1) = 1;
%     x(x>sx) = sx;
%     y(y<1) = 1;
%     y(y>sy) = sy;
%     z(z<1) = 1;
%     z(z>sz) = sz;



%     figure(1)
%     subplot(5,1,1)
%     imagesc([I2(:,:,round(sz/2)),I0(:,:,round(sz/2)),dI(:,:,round(sz/2))])
%     colormap gray
%     axis image
%     brighten(0.4)
%     drawnow
%     subplot(5,1,2)
%     plot(E)

   
end
end

function input = fft3(input)
%imspace = fftshift(imspace,3);
input = fft(input,[],1);
input = fft(input,[],2);
input = fft(input,[],3);
%input = fft(fft(fft(input,[],1),[],2),[],3);
end

function input = ifft3(input)

%input = ifft(ifft(ifft(input,[],1),[],2),[],3);
input = ifft(input,[],1);
input = ifft(input,[],2);
input = ifft(input,[],3);
%imspace = fftshift(imspace,3);

end

function input = fftshift3(input)

input = fftshift(input,1);
input = fftshift(input,2);
input = fftshift(input,3);
end

function input = ifftshift3(input)

input = ifftshift(input,1);
input = ifftshift(input,2);
input = ifftshift(input,3);
end