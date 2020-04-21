function llr = Patch_tracking_3D_with_guide_random_search(Image,Nphase,offset_init,patch_size,search_size,patch_shift)

[sx,sy,nof,nslice] = size(Image);
% Nphase = 5;%size(para.Recon.bins,1);
Ncycle = nof/Nphase;
im_bin = reshape(Image,[sx,sy,Nphase,Ncycle,nslice]);
im_bin = permute(im_bin,[1,2,5,3,4]);

if ~exist('patch_size')
    patch_size = [5,5,2];
end
if ~exist('search_size')
    search_size = [9,9,4];
end
if ~exist('patch_shift')
    patch_shift = [5,5,2];
end

x_begin_all = 1:patch_shift(1):sx-patch_size(1)+1;
y_begin_all = 1:patch_shift(2):sy-patch_size(2)+1;
z_begin_all = 1:patch_shift(3):nslice-patch_size(3)+1;
% c_begin_all = 1:Nphase:Nphase-patch_size(4)+1;

if z_begin_all(end) + patch_size(3)-1 < nslice
    z_begin_all(end+1) = nslice-patch_size(3) + 1;
end
if x_begin_all(end) + patch_size(1)-1 < sx
    x_begin_all(end+1) = sx-patch_size(1) + 1;
end
if y_begin_all(end) + patch_size(2)-1 < sy
    y_begin_all(end+1) = sy-patch_size(2) + 1;
end
Nx = length(x_begin_all);
Ny = length(y_begin_all);
Nz = length(z_begin_all);
% Nc = length(c_begin_all);

%Npatch_all = Nx*Ny*Nz*Nc;
%int_patch = zeros([patch_size,Npatch_all]);
std_map = std(im_bin,1,5);
int_map = sum(abs(im_bin),5);
%std_patch = zeros([patch_size,Npatch_all]);

%offset_init = squeeze(sum(offset_init,2));
offset_init = reshape(offset_init,[3,Nphase,Ncycle]);
offset_init = offset_init - offset_init(:,:,1);
offset_init(:,:,2:end) = diff(offset_init,1,3);
% offset_init = squeeze(mean(offset_init,2));
% offset_init(3,:) = offset_init(3,:)/3;
offset_init = round(offset_init);
% offset_init = cat(1,offset_init,zeros(1,Ncycle));
% offset_init = permute(offset_init,[2,1]);
offset_init = -offset_init;



offset_all = [];
for iphase = 1:Nphase
    clear idx_all int_patch std_patch offset_all_phase
    N = 0;
    Image_temp = squeeze(im_bin(:,:,:,iphase,:));
    offset_phase = squeeze(offset_init(:,iphase,:))';
    std_temp = std_map(:,:,:,iphase);
    int_temp = int_map(:,:,:,iphase);
    
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                x_temp = x_begin_all(i):x_begin_all(i)+patch_size(1)-1;
                y_temp = y_begin_all(j):y_begin_all(j)+patch_size(2)-1;
                z_temp = z_begin_all(k):z_begin_all(k)+patch_size(3)-1;
                %                 c_temp = c_begin_all(l):c_begin_all(l)+patch_size(4)-1;
                [x_temp,y_temp,z_temp] = ndgrid(x_temp,y_temp,z_temp);
                idx_temp = sub2ind([sx,sy,nslice],x_temp,y_temp,z_temp);
                
                N = N+1;
                %idx_all(N) = sum(vec(mask(x_temp,y_temp,z_temp)));
                idx_all(:,:,:,N) = idx_temp;
                int_patch(:,:,:,N) = int_temp(idx_temp);
                std_patch(:,:,:,N) = std_temp(idx_temp);
            end
        end
    end
    
    
    idx = sum(sum(sum(abs(int_patch))));
    idx_std = sum(sum(sum(abs(std_patch))));
    keep = (idx_std > max(idx_std)/10) & (idx > max(idx)/10);
    
    idx_all = idx_all(:,:,:,keep);
    idx_add = (0:Ncycle-1)*sx*sy*nslice;
    idx_all = idx_all + permute(idx_add,[1,3,4,5,2]);
    
    Npatch = sum(keep(:));
    idx_all = reshape(idx_all,[patch_size,Npatch*Ncycle]);
    
    [~,searching_order] = sort(rand([1,Npatch*Ncycle]));
    %keep = idx_all~=0;
    
    %% patch search
    
    mask = zeros(size(Image_temp));
    n=0;
    tic
    for i=1:Npatch*Ncycle
        %     if mod(i,Npatch)==0
        %         fprintf([num2str(i/(Npatch*Ncycle)),'\n'])
        %     end
        order_temp = searching_order(i);
        idx_temp = idx_all(:,:,:,order_temp);
        flag = sum(vec(mask(idx_temp) == 0));
        
        if flag
            n = n+1;
            [x_temp,y_temp,z_temp,c_temp] = ind2sub([sx,sy,nslice,Ncycle],idx_temp);
            x_begin = x_temp(1);
            y_begin = y_temp(1);
            z_begin = z_temp(1);
            c_temp = c_temp(1);
            mask(idx_temp) = mask(idx_temp) + 1;
            
            offset_all_temp = zeros(3,Ncycle);
            offset_all_temp(:,c_temp) = [x_begin,y_begin,z_begin];
            if c_temp~=1
                for icycle = c_temp-1:-1:1
                    offset_temp = -offset_phase(icycle+1,:);
                    
                    x_window_begin = max(x_begin-round((search_size(1)-patch_size(1))/2)+offset_temp(1),1); x_window_end = min(x_window_begin-1+search_size(1),sx);
                    if x_window_begin==1
                        x_window_end = max(x_window_end,x_window_begin+patch_size(1)-1);
                    end
                    if x_window_end==sx
                        x_window_begin = min(x_window_begin,x_window_end-patch_size(1)+1);
                    end
                    
                    y_window_begin = max(y_begin-round((search_size(2)-patch_size(2))/2)+offset_temp(2),1); y_window_end = min(y_window_begin-1+search_size(2),sy);
                    if y_window_begin==1
                        y_window_end = max(y_window_end,y_window_begin+patch_size(2)-1);
                    end
                    if y_window_end==sy
                        y_window_begin = min(y_window_begin,y_window_end-patch_size(2)+1);
                    end
                    
                    z_window_begin = max(z_begin-round((search_size(3)-patch_size(3))/2)+offset_temp(3),1); z_window_end = min(z_window_begin-1+search_size(3),nslice);
                    if z_window_begin==1
                        z_window_end = max(z_window_end,z_window_begin+patch_size(3)-1);
                    end
                    if z_window_end==nslice
                        z_window_begin = min(z_window_begin,z_window_end-patch_size(3)+1);
                    end
                    
                    patch_temp = Image_temp(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle+1);
                    search_window = Image_temp(x_window_begin:x_window_end,y_window_begin:y_window_end,z_window_begin:z_window_end,icycle);
                    
                    [offset] = patch_search_one_3D(patch_temp,search_window);
                    
                    begin_all(1) = x_window_begin-1+offset(1); x_begin = begin_all(1);
                    begin_all(2) = y_window_begin-1+offset(2); y_begin = begin_all(2);
                    begin_all(3) = z_window_begin-1+offset(3); z_begin = begin_all(3);
                    
                    offset_all_temp(:,icycle) = begin_all;
                    mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle) = mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle) + 1;
                end
            end
            if c_temp~=Ncycle
                x_begin = x_temp(1);
                y_begin = y_temp(1);
                z_begin = z_temp(1);
                
                for icycle = c_temp:1:Ncycle-1
                    offset_temp = offset_phase(icycle+1,:);
                    
                    x_window_begin = max(x_begin-1+offset_temp(1),1); x_window_end = min(x_begin-1+patch_size(1)+1+offset_temp(1),sx);
                    if x_window_begin==1
                        x_window_end = max(x_window_end,x_window_begin+patch_size(1)-1);
                    end
                    if x_window_end==sx
                        x_window_begin = min(x_window_begin,x_window_end-patch_size(1)+1);
                    end
                    
                    y_window_begin = max(y_begin-1+offset_temp(2),1); y_window_end = min(y_begin-1+patch_size(2)+1+offset_temp(2),sy);
                    if y_window_begin==1
                        y_window_end = max(y_window_end,y_window_begin+patch_size(2)-1);
                    end
                    if y_window_end==sy
                        y_window_begin = min(y_window_begin,y_window_end-patch_size(2)+1);
                    end
                    
                    z_window_begin = max(z_begin-1+offset_temp(3),1); z_window_end = min(z_begin-1+patch_size(3)+1+offset_temp(3),nslice);
                    if z_window_begin==1
                        z_window_end = max(z_window_end,z_window_begin+patch_size(3)-1);
                    end
                    if z_window_end==nslice
                        z_window_begin = min(z_window_begin,z_window_end-patch_size(3)+1);
                    end
                    
                    patch_temp = Image_temp(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle);
                    search_window = Image_temp(x_window_begin:x_window_end,y_window_begin:y_window_end,z_window_begin:z_window_end,icycle+1);
                    
                    [offset] = patch_search_one_3D(patch_temp,search_window);
                    
                    begin_all(1) = x_window_begin-1+offset(1); x_begin = begin_all(1);
                    begin_all(2) = y_window_begin-1+offset(2); y_begin = begin_all(2);
                    begin_all(3) = z_window_begin-1+offset(3); z_begin = begin_all(3);
                    
                    offset_all_temp(:,icycle+1) = begin_all;
                    mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle+1) = mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,z_begin:z_begin+patch_size(3)-1,icycle+1) + 1;
                end
            end
            offset_all_phase(:,:,n) = offset_all_temp;
        end
    end
    
    offset_all = cat(3,offset_all,offset_all_phase);
    phase_all(iphase) = n;
    mask_all(:,:,:,iphase,:) = mask;
    toc
end

n = sum(phase_all);

phase_idx = [];
for i=1:Nphase
    phase_idx = cat(2,phase_idx,ones(1,phase_all(i))*i);
end

patch_begin_all = offset_all;
patch_end_all = offset_all + reshape(patch_size(1:3),3,1) - 1;

idx = zeros([patch_size,Ncycle,n]);
for i=1:n
    %     fprintf([num2str(i/n),'\n'])
    for j=1:Ncycle
        x_temp = patch_begin_all(1,j,i):patch_end_all(1,j,i);
        y_temp = patch_begin_all(2,j,i):patch_end_all(2,j,i);
        z_temp = patch_begin_all(3,j,i):patch_end_all(3,j,i);
        p_temp = phase_idx(i);
        [x_temp,y_temp,z_temp,p_temp] = ndgrid(x_temp,y_temp,z_temp,p_temp);
        idx_temp = sub2ind([sx,sy,nslice,Nphase],x_temp,y_temp,z_temp,p_temp);
        idx(:,:,:,j,i) = idx_temp;
    end
end
idx = reshape(idx,prod(patch_size),Ncycle,n);
idx_add = sx*sy*nslice*Nphase*(0:Ncycle-1);
idx = idx+idx_add;

mask = reshape(mask_all,[sx,sy,nslice,Nphase*Ncycle]);
mask_LR = mask;
mask_LR(mask_LR==0) = 1;
mask_LR = 1./mask_LR;
mask = logical(mask);
Npatch = n;

llr.idx = idx;
llr.Npatch = Npatch;
llr.mask_intensity = mask_LR;
llr.mask = mask;


return










