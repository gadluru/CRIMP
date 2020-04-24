function step_size = line_search(old, update, Data, para)
%--------------------------------------------------------------------------
%   [step] = line_search(old, update, Data, para)
%--------------------------------------------------------------------------
%   Line search called in a conjugate gradient algorithm
%--------------------------------------------------------------------------
%   Inputs:      
%       - old       [sx, sy, nof, ...]
%       - update    [sx, sy, nof, ...]
%       - Data      [structure]
%       - para      [structure]
%               
%       - old       image from previous iteration
%       - update    update term
%       - Data      see 'help STCR_conjugate_gradient.m'
%       - para      see 'help STCR_conjugate_gradient.m'
%--------------------------------------------------------------------------
%   Output:
%       - step      [scalar]
%
%       - step      step size for CG update
%--------------------------------------------------------------------------
%   This function trys to find a suitable step size to perform a CG update.
%   The function starts with a step size adopted from last iteration, and
%   multiply it by 1.3 (magic number). If the step size yeilds a cost that
%   is larger than the previous cost, it shrinks the step size by 0.8
%   (magic number again). If it yeilds a cost that is smaller than the
%   previous cost, it will increase the step size by 1.3 until it no longer
%   yeild a smaller cost. The maximum number of trys is 15.
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------
step_start = para.Recon.step_size(end)*1.3;%magic number
tau = 0.8;
max_try = 15;
step_size = step_start;

cost_old = para.Cost.totalCost(end);
flag = 0;

for i=1:max_try
%     fprintf(['Iter = ' num2str(i) '...\n'])
    
    new = old + step_size*update;
    fidelity_new = compute_fidelity_for_line_search_yt(new,Data,para);
    cost_new = Cost_STCR(fidelity_new,new,para.Recon.weight_sTV,para.Recon.weight_tTV);
    
    if cost_new > cost_old && flag == 0
        step_size = step_size*tau;
    elseif cost_new < cost_old 
        step_size = step_size*1.3;
        cost_old = cost_new;
        flag = 1;
    elseif cost_new > cost_old && flag == 1
        step_size = step_size/1.3;
%         fprintf(['Step = ' num2str(step) '...\n'])
%         fprintf(['Cost = ' num2str(round(cost_old)) '...\n'])
        return
    end
end
% fprintf(['Step = ' num2str(step) '...\n'])
% fprintf(['Cost = ' num2str(round(cost_new)) '...\n'])
end
