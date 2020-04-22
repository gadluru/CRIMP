function update = patch_low_rank_only(Image,llr)

Image = permute(Image,[1,2,4,3]);
patch_all = Image(llr.idx);
update = zeros(size(Image),'like',Image);

for i=1:llr.Npatch
    [u,s,v] = svd(patch_all(:,:,i),0);
    s = s-20;
    s(s<0) = 0;
    update(llr.idx(:,:,i)) = update(llr.idx(:,:,i)) + u*s*v';
end

update = update.*llr.mask_intensity; 
update = permute(update,[1,2,4,3]);