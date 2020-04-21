function [Data,para] = get_Data_SMS(Data,para)

Data.kSpace = fftshift2(Data.kSpace);
Data.mask = logical(abs(Data.kSpace(:,:,:,1,:,:,:)));
Data.kSpace = ifft2(Data.kSpace);
Data.kSpace = fftshift2(Data.kSpace);
Data.kSpace = fft2(Data.kSpace);
Data.kSpace = Data.kSpace.*Data.mask;
nSMS = para.Recon.nSMS;
Data.SMS = exp(1i*(0:nSMS-1)*2*pi/nSMS.*(0:nSMS-1).');
Data.SMS = permute(Data.SMS,[3,4,5,6,1,7,2]);
kSpace_sms = sum(Data.kSpace.*conj(Data.SMS),7);
Data.first_est = ifft2(kSpace_sms);
if ~isfield(Data,'sens_map')
    Data.sens_map = get_sens_map(Data.first_est,'SMS');
end
para.Recon.kSpace_size = [size(Data.kSpace,1),size(Data.kSpace,2)];
para.Recon.image_size = [size(Data.kSpace,1),size(Data.kSpace,2)];
para.Recon.sx = size(Data.kSpace,1);
Data.filter = ramp_filter_for_pre_interp(para);
Data.first_est = ifft2(kSpace_sms.*Data.filter);
Data.first_est = Data.first_est.*Data.sens_map;
Data.first_est = sum(Data.first_est,4);
