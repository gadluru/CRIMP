function update = patch_ttv(Image,llr,para)

Image = permute(Image,[1,2,4,3]);
patch_all = Image(llr.idx);

update_tv = compute_tTV_yt(permute(patch_all,[1,3,2]),para.Recon.weight_tTV,1e-7);
update_tv = permute(update_tv,[1,3,2]);
update_t = zeros(size(Image),'like',Image);

for i=1:llr.Npatch
    update_t(llr.idx(:,:,i)) = update_t(llr.idx(:,:,i)) + update_tv(:,:,i);
end

update = update_t.*llr.mask_intensity;    
update = permute(update,[1,2,4,3]);