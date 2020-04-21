function fidelity_norm = compute_fidelity_for_line_search_yt(image,Data,para)

switch para.Recon.type
        
    case 'seperate SMS test'
        fidelity_norm = 0;
        for i=1:para.Recon.no_comp
            fidelity_update = image.*Data.sens_map(:,:,:,i,:,:);
            fidelity_update = fft2(fidelity_update);
            fidelity_update = sum(fidelity_update.*Data.SMS,5);
            fidelity_update = fidelity_update(Data.mask);
            fidelity_update = Data.kSpace(:,i) - fidelity_update;
            fidelity_norm = fidelity_norm + sum(abs(fidelity_update(:)).^2/prod(para.Recon.kSpace_size));
        end
        fidelity_norm = sqrt(fidelity_norm);
        return
        
end
        
end