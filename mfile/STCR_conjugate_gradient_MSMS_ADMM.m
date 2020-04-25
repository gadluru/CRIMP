function [Image,para] = STCR_conjugate_gradient_MSMS_ADMM(Data,para)
%--------------------------------------------------------------------------
%   [Image,para] = STCR_conjugate_gradient_MSMS_ADMM(Data,para)
%--------------------------------------------------------------------------
%   Solve MRI reconstruction problem within the ADMM iteration using
%   conjugated gradient algorithm.
%--------------------------------------------------------------------------
%   Inputs:
%       - Data                      [structure] 
%           Data.kSpace             [nm, nc]
%           Data.sens_map           [1,  1,   1,   ns, nc]
%           Data.first_est          [sx, sy,  nof, ns]
%           Data.first_guess        [sx, sy,  nof, ns]
%           Data.mask
%
%               'sx'    number of readout point along a ray
%               'sy'    for radial k-space, same as sx
%               'nm'    number of total non-zero measurements
%               'nor'   number of rays per time frame
%               'nof'   number of time frames
%               'nc'    number of coils
%           
%       - para                      [structure]
%           para.setting            [structure]
%               setting.ifplot      [0 or 1]
%               setting.ifGPU       [0 or 1]
%           para.Recon              [structure]
%               Recon.weight_tTV    [scalar]
%               Recon.weight_sTV    [scalar]
%               Recon.weight_l2     [scalar]
%               Recon.epsilon       [scalar]         
%               Recon.step_size     [scalar]
%
%       - Data
%           Data.kSpace             measured k-space data "d"
%           Data.sens_map           sensitivity map
%           Data.first_est          initial estimation of "x"
%           Data.first_guess        "Phat"
%
%       -para
%           para.setting.ifplot     display reconstruction process
%           para.setting.ifGPU      run function on a NVIDIA GPU
%           para.Recon.weight_tTV   "lambda_t"
%           para.Recon.weight_sTV   "lambda_s"
%           para.Recon.weight_l2    "rho"
%           para.Recon.epsilon      "epsilon"
%           para.Recon.step_size    initial CG update step size
%--------------------------------------------------------------------------
%   Output:
%       - Image     [sx, sy, nof, ...]
%       - para      [structure]
%
%       - Image     reconstructed images "m"
%--------------------------------------------------------------------------
%   A standard cost function it solves is the spatially and temporally
%   constrained reconstruction (STCR):
%
%   || Am - d ||_2^2 + lambda_t || TV_t m ||_1 
%                    + lambda_t sum_i || P_i m ||_1
%                    + rho || Phat - m - Y ||_2^2
%
%   "A"         sampling matrix includes sensitivity maps, Fourier 
%               transform, and undersampling mask
%   "m"         image to be reconstructed
%   "d"         measured k-space data
%   ||.||_2^2   l2 norm
%   ||.||_1     l1 norm
%   "lambda_t"  temporal constraint weight
%   TV_t        temporal total variation (TV) operator (finite difference)
%               sqrt( abs(m_t+1 - m_t)^2 + epsilon )
%   "epsilon"   small term to aviod singularity
%   "rho"       l2 norm constraint weight
%   sum_i       sum across all patches
%   P_i         patch extraction operator
%   "Phat"      low-rank images
%   "Y"         ADMM Lagrangian multiplier
%--------------------------------------------------------------------------
%   Reference:
%       [1] Whole-heart, ungated, free-breathing, cardiac-phase-resolved 
%           myocardial perfusion MRI by using Continuous Radial Interleaved
%           simultaneous Multi-slice acquisitions at sPoiled steady-state 
%           (CRIMP). MRM, in press.
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------
disp('Performing iterative CG reconstruction...');
disp('Showing progress...')

ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
epsilon        = para.Recon.epsilon;
weight_l2      = para.Recon.weight_l2;
para.Recon.step_size = para.Recon.step_size(1);


if isfield(Data,'first_est')
    new_img_x = single(Data.first_est);
else
    new_img_x = Data.first_guess;
end

if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if ifGPU
%    Data.kSpace        = gpuArray(Data.kSpace);
    new_img_x          = gpuArray(new_img_x);
    Data.sens_map      = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    if isfield(Data,'mask')
        Data.mask          = gpuArray(Data.mask);
    end
    if isfield(Data,'filter')
        Data.filter        = gpuArray(Data.filter);
    end
%     epsilon = gpuArray(epsilon);
end

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'l2Norm',[],'totalCost',[]);

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
temporal = @(im) compute_tTV_yt(im,weight_tTV,epsilon);

for iter_no = 1:para.Recon.noi
    t1 = tic;

%% fidelity term/temporal/spatial TV
    [update_term,fidelity_norm] = fidelity(sr_forward(new_img_x,Data.size));
    update_term = sr_backward(update_term,Data.size);
    update_term = update_term + temporal(new_img_x)*0.5;

    update_patch = patch_ttv(gather(new_img_x),Data.llr,para)*0.5;
    update_bins = compute_tTV_bins_no_lr(gather(new_img_x),weight_tTV,epsilon,para.Recon.bins)*0.5;
    update_bins(permute(Data.llr.mask,[1,2,4,3])) = update_patch(permute(Data.llr.mask,[1,2,4,3]));
    update_term = update_term + update_bins; clear update_bins update_patch
    if isfield(Data, 'Y')
        update_term = update_term + (Data.first_guess - new_img_x - Data.Y)*weight_l2;
    end

%% conjugate gradient
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%% line search
    if isfield(Data, 'Y')
        para.Cost = Cost_STCR_step_3(fidelity_norm, gather(new_img_x), weight_sTV, weight_tTV, Data.llr, Data.first_guess - Data.Y, weight_l2, para.Cost);
    else
        para.Cost = Cost_STCR_step_2(fidelity_norm, new_img_x, weight_sTV, weight_tTV, Data.llr, para.Cost);
    end
    clear fidelity_update
    step_size = line_search_step_2(new_img_x,update_term_old,Data,para);
    para.Recon.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size * update_term_old;
    update_term_old = gather(update_term_old);

%% show reconstruction progress
    if ifplot == 1
        showImage3D(new_img_x,para.Cost)
    end

%% stoping creteria
    if para.Recon.break && iter_no > 1
        if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end
    
    fprintf(['Iteration = ' num2str(iter_no) '...']);
    para.Recon.time(iter_no) = toc(t1);
    toc(t1)
end

Image = squeeze(gather(new_img_x));
para.Recon.time_total = sum(para.Recon.time);
figure, plotCost(para.Cost); drawnow
fprintf(['Iterative reconstruction running time is ' num2str(para.Recon.time_total) 's' '\n'])
end