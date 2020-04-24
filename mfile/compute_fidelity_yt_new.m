function [fidelity_update, fidelity_norm] = compute_fidelity_yt_new(image, Data, para)
%--------------------------------------------------------------------------
%   [fidelity_update, fidelity_norm] 
%       = compute_fidelity_yt_new(image, Data, para)
%--------------------------------------------------------------------------
%   Compute fidelity update term of a MRI reconstruction problem
%--------------------------------------------------------------------------
%   Inputs:      
%       - image             [sx, sy, nof, ...]
%       - Data              [structure]
%       - para              [structure]
%
%       - image             image
%       - Data              see 'help STCR_conjugate_gradient.m'
%       - para              see 'help STCR_conjugate_gradient.m'
%--------------------------------------------------------------------------
%   Output:
%       - fidelity_update   [sx, sy, nof, ...]
%       - fidelity_norm     [scalar]
%
%       - fidelity_update   A^H (Am - d)
%       - fidelity_norm     || Am - d ||_2^2
%--------------------------------------------------------------------------
%   A standard fidelity term it solves is:
%
%   || Am - d ||_2^2
%
%   and the output is:
%
%   A^H (Am - d)
%
%   see 'help STCR_conjugate_gradient.m' for more information.
%--------------------------------------------------------------------------
%   Reference:
%       [1]     Acquisition and reconstruction of undersampled radial data 
%               for myocardial perfusion MRI. JMRI, 2009, 29(2):466-473.
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------
switch para.Recon.type
   
    case 'seperate SMS test'

        fidelity_update_all = single(zeros(size(image),class(image)));
        fidelity_norm = 0;
        for i=1:para.Recon.no_comp
            fidelity_update = image.*Data.sens_map(:,:,:,i,:,:);
            fidelity_update = fft2(fidelity_update);
            fidelity_update = sum(fidelity_update.*Data.SMS,5);
            siz = size(fidelity_update);
            fidelity_update_temp = zeros(siz,class(fidelity_update));
            fidelity_update = fidelity_update(Data.mask);
            fidelity_update = Data.kSpace(:,i) - fidelity_update;
            fidelity_norm = fidelity_norm + sum(abs(fidelity_update(:)).^2/prod(para.Recon.kSpace_size));
            fidelity_update_temp(Data.mask) = fidelity_update;

            fidelity_update = fidelity_update_temp;clear fidelity_update_temp
            fidelity_update = sum(fidelity_update.*conj(Data.SMS),7);
            if isfield(Data,'filter')
                fidelity_update = bsxfun(@times,fidelity_update,Data.filter);
            end
            int_mask = sum(Data.mask,7);
            int_mask(int_mask==0) = 1;
            int_mask = 1./int_mask;
            fidelity_update = fidelity_update.*int_mask;
            
            fidelity_update = ifft2(fidelity_update);
            fidelity_update_all = fidelity_update_all + bsxfun(@times,fidelity_update,Data.sens_map_conj(:,:,:,i,:,:));
        end
        fidelity_update = fidelity_update_all;
        fidelity_norm = sqrt(fidelity_norm);
        return
        
        
  
end
