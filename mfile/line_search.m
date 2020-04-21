function step_size = line_search(old,update,Data,para)

step_start = para.step_size(end)*1.3;%magic number
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
