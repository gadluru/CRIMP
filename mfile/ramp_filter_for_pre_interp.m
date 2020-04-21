function filter = ramp_filter_for_pre_interp(para)

[X,Y] = meshgrid(1:para.Recon.kSpace_size(1),1:para.Recon.kSpace_size(2));
X = X - (para.Recon.kSpace_size(1)+1)/2;
Y = Y - (para.Recon.kSpace_size(2)+1)/2;
filter = sqrt(X.^2 + Y.^2);
fully_sampled_radius = para.Recon.nor*2*para.over_sampling/pi;
filter(filter<fully_sampled_radius) = fully_sampled_radius;
filter = filter/(para.Recon.sx+1)*2;
filter(filter>1) = 1;
filter = fftshift2(filter);
filter = single(filter);

end