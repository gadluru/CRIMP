ccc
temp = dir('./RawData/meas_MID00051_FID30521_UCAIR_Ungated2D_SPGR_MB3_11_3ml_stress_Kspace_PCA.mat');
if isempty(temp)
    temp = dir('./RawData/kSpace_1.mat');
    if ~isempty(temp)
        load('./RawData/kSpace_1.mat')
    else
        error('missing file kSpace_1.mat under ./RawData')
    end
    temp = dir('./RawData/kSpace_2.mat');
    if ~isempty(temp)
        load('./RawData/kSpace_2.mat')
    else
        error('missing file kSpace_2.mat under ./RawData')
    end
    temp = dir('./RawData/kSpace_3.mat');
    if ~isempty(temp)
        load('./RawData/kSpace_3.mat')
    else
        error('missing file kSpace_3.mat under ./RawData')
    end
    temp = dir('./RawData/kSpace_4.mat');
    if ~isempty(temp)
        load('./RawData/kSpace_4.mat')
    else
        error('missing file kSpace_4.mat under ./RawData')
    end
    temp = dir('./RawData/kSpace_info.mat');
    if ~isempty(temp)
        load('./RawData/kSpace_info.mat')
    else
        error('missing file kSpace_info.mat under ./RawData')
    end
    kSpace = cat(2, kSpace_1, kSpace_2, kSpace_3, kSpace_4);
    save('./RawData/meas_MID00051_FID30521_UCAIR_Ungated2D_SPGR_MB3_11_3ml_stress_Kspace_PCA.mat', 'kSpace', 'kSpace_info')
    !rm ./RawData/kSpace*
end