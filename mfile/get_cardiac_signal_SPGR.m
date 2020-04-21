function [cardiac_signal,resp_signal,mask] = get_cardiac_signal_SPGR(Image)
% [cardiac_signal,resp_signal,mask] = get_cardiac_signal_SPGR(Image)
% 
% get cardiac and respiratory signals from reference images

%% cardiac signal
image_for_std = crop_half_FOV(abs(Image(:,:,1:28,:)));
nslice = size(image_for_std,4);
sx = size(image_for_std,1);
nof = size(Image,3);
image_for_std = reshape(image_for_std,[sx,sx,7,4,nslice]);
std_map = std(image_for_std,0,3);
std_map = squeeze(sum(std_map,4));
std_map = imgaussfilt(std_map,5);
filter = fspecial('gaussian',sx,sx/10);
std_map = std_map.*filter;
mask = zeros(sx,sx,nslice);
for i=1:nslice
    [x,y] = find(std_map(:,:,i)==max(max(std_map(:,:,i))));
    mask(x,y,i) = 1;
    mask(:,:,i) = bwdist(mask(:,:,i));
end

mask = mask<45;


dt = 35;
Fs = 1000/dt;
df = Fs/nof;

fpass_cardiac = [0.5/df/(nof/2),2.2/df/(nof/2)];


cardiac_signal = permute(mask,[1,2,4,3]).*crop_half_FOV(abs(Image(:,:,:,:)));
cardiac_signal = permute(cardiac_signal,[1,2,4,3]);
cardiac_signal = reshape(cardiac_signal,[sx*sx*nslice,nof]);
cardiac_signal = sum(cardiac_signal,1);
cardiac_signal = bandpass(cardiac_signal,fpass_cardiac);
%% respiration signal 

resp_range = round((0.2/df):(0.5/df));
fpass_resp = [0.2/df/(nof/2),0.5/df/(nof/2)];

if resp_range==0
    resp_signal = ones(1,nof);
    return
end

resp_signal = squeeze(sum(abs(Image(:,:,:,:)),2));
resp_signal = permute(resp_signal,[1,3,2]);
resp_signal = reshape(resp_signal,[sx*2*nslice,nof]);
resp_coeff = pca(resp_signal);
resp_signal_fft = fft(resp_coeff(:,1:10),[],1);

[~,idx] = max(max(abs(resp_signal_fft(resp_range,:))),[],2);

[~,idx_max] = max(abs(resp_signal_fft(resp_range,idx)),[],1);
peak_frequency = (resp_range(1) + idx_max - 1)*df/Fs;

resp_signal = lowpass(resp_coeff(:,idx),peak_frequency);

end