function Image = LLR_patch_tracking_4D_ADMM(Data,para)

disp('Performing ADMM reconstruction...');
disp('Showing progress...')

para.Recon.weight_l2 = 0.01;
Image = STCR_conjugate_gradient_MSMS_ADMM(Data,para);

Data.Y = zeros(size(Data.first_guess));

for ADMM_iter = 1:2
    %% ADMM update step 2
    % motion tracked patches
    Image_mc_lr = patch_low_rank_only(Image,Data.llr);
    Image_mc_lr = permute(Image_mc_lr,[1,2,4,3]);
    
    % not moving patches
    Image_lr = reshape(Image,[para.Recon.sx,para.Recon.sx,Data.size.Nphase,Data.size.Ncycle,Data.size.nset*Data.size.nSMS]);
    Image_lr = permute(Image_lr,[1,2,5,4,3]);
    Image_lr = LowRank_patch_yt(Image_lr);
    Image_lr = permute(Image_lr,[1,2,5,4,3]);
    Image_lr = reshape(Image_lr,[para.Recon.sx,para.Recon.sx,Data.size.Nphase*Data.size.Ncycle,Data.size.nset*Data.size.nSMS]);
    Image_lr = permute(Image_lr,[1,2,4,3]);
    
    % combine two LR images
    Image_lr(Data.llr.mask) = Image_mc_lr(Data.llr.mask);
    Image_lr = permute(Image_lr,[1,2,4,3]);
    
    %% ADMM update step 3
    Data.first_guess = Image_lr;
    Data.Y = Data.Y + Image - Image_lr;
    
    %% ADMM update step 1
    % CG solve for constrained problem:
    % min_x = ||Ax - y||_2^2 + lt||dt x||_1 + lt||dt p(x)||_1 + l2||x - x^||_2
    Image = STCR_conjugate_gradient_MSMS_ADMM(Data,para);
end


end

