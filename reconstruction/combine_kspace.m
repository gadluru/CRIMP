ccc
load('./RawData/kSpace_1.mat')
load('./RawData/kSpace_2.mat')
load('./RawData/kSpace_3.mat')
load('./RawData/kSpace_info.mat')
kSpace = cat(2, kSpace_1, kSpace_2, kSpace_3);
save('./RawData/meas_MID00051_FID30521_UCAIR_Ungated2D_SPGR_MB3_11_3ml_stress_Kspace_PCA.mat', 'kSpace', 'kSpace_info')
!rm ./RawData/kSpace*