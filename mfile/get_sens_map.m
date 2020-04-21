function sens_map = get_sens_map(im,options)

smooth = 10;

if size(im,4) == 1 & ~contains(options,'3D')
    sens_map = ones(size(im,1),size(im,2),'single');
    return
end
    
switch options
    case '2D'
        im_for_sens = squeeze(sum(im,3));
        sens_map(:,:,1,:) = ismrm_estimate_csm_walsh_optimized_yt(im_for_sens,smooth);
    case '3D'
        [sx,sy,sz,~,coils] = size(im);
        im_for_sens = squeeze(sum(im,4));
        sens_map = single(zeros(sx,sy,sz,1,coils));
        for i=1:sz
            sens_map(:,:,i,1,:) = ismrm_estimate_csm_walsh_optimized_yt(squeeze(im_for_sens(:,:,i,:)),smooth);
        end
    case 'SMS'
        [sx,sy,nof,coils,nSMS,ns] = size(im);
        im = reshape(im,[sx,sy,nof,coils,nSMS*ns]);
        im_for_sens = squeeze(sum(im,3));
        sens_map = zeros(sx,sy,1,coils,nSMS*ns);
        for i=1:nSMS*ns
            sens_map(:,:,1,:,i) = ismrm_estimate_csm_walsh_optimized_yt(im_for_sens(:,:,:,i),smooth);
        end
        sens_map = reshape(sens_map,[sx,sy,1,coils,nSMS,ns]);
    case 'ESPIRIT'
        run('/v/raid1b/ytian/From_raid1a/mfile_from_else/ESPIRiT/setPath.m')
        eigThresh_k = 0.02; % threshold of eigenvectors in k-space
        eigThresh_im = 0.9; % threshold of eigenvectors in image space
        kernel_size = [6,6];
        center_rad = 10;
        [sx,sy,~,coils,slice] = size(im);
        for i=1:slice
            im_for_sens = squeeze(mean(im(:,:,:,:,i),3));
            im_for_sens = fftshift(fftshift(im_for_sens,1),2);
            k_for_sens = fft2(im_for_sens);
            k_for_sens = fftshift(fftshift(k_for_sens,1),2);
            k_center = round(size(k_for_sens,1)/2);
            k_for_sens = k_for_sens(k_center-center_rad+1:k_center+center_rad,k_center-center_rad+1:k_center+center_rad,:);
            
            [k,S] = dat2Kernel(k_for_sens,kernel_size);
            idx = find(S >= S(1)*eigThresh_k,1,'last');
            [M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);
            
            weights = (W - eigThresh_im)./(1-eigThresh_im).* (W> eigThresh_im);
            weights = -cos(pi*weights)/2 + 1/2;
            weights(weights==0) = 0.01;
            sens_map_temp = weights.*permute(M,[1,2,4,3]);
            sens_map_temp = sens_map_temp(:,:,8,:);
            sens_map(:,:,:,:,i) = permute(sens_map_temp,[1 2 5 4 3]);
        end
        %sens_map = sum(sens_map,3);
        
        
end

sens_map_scale = max(abs(sens_map(:)));
sens_map = sens_map/sens_map_scale;
sens_map_conj = conj(sens_map);
switch options
    case '3D'
        sens_correct_term = 1./sum(sens_map_conj.*sens_map,5);
    case 'ESPIRIT'
        sens_correct_term = 1./sum(sum(sens_map_conj.*sens_map,4),5);
        sens_correct_term = 1;
    otherwise
        sens_correct_term = 1./sum(sens_map_conj.*sens_map,4);
end
sens_correct_term = sqrt(sens_correct_term);
sens_map = bsxfun(@times,sens_correct_term,sens_map);

sens_map = single(sens_map);