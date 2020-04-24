function sens_map = get_sens_map(im, options)

smooth = 10;

if size(im,4) == 1 && ~contains(options,'3D')
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
        
end

sens_map_scale = max(abs(sens_map(:)));
sens_map = sens_map/sens_map_scale;
sens_map_conj = conj(sens_map);
switch options
    case '3D'
        sens_correct_term = 1./sum(sens_map_conj.*sens_map,5);
    otherwise
        sens_correct_term = 1./sum(sens_map_conj.*sens_map,4);
end
sens_correct_term = sqrt(sens_correct_term);
sens_map = bsxfun(@times,sens_correct_term,sens_map);

sens_map = single(sens_map);