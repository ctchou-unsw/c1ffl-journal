% To plot the graph for intermediate approximation for the paper 

% Specify test number 
test_num = 1;

% Load data 
eval(['load eval_int_approx_v1_data_ex',int2str(test_num)])

% Plot: Style 1: With only one L(t) and one Lhat(t) 
if test_num == 1
    % A figure with only two curves  
    sim_index = 1;
    figure(1)
    plot(vec_time,mat_llr_exact(:,sim_index),'b-', ...
         vec_time,mat_llr_appro(:,sim_index),'r-', ...
         'Linewidth',2);  
    legend({'L(t)','$\hat{L}(t)$'},'Location','NorthWest', ...
            'FontWeight','Bold','FontSize',20,'Interpreter','latex')
    xlabel('time','FontSize',20) 
    ylabel('Log-likelihood Ratio','FontSize',20)
    eval(['print -depsc plot_ia_2curves_',int2str(test_num)])
end 

% Plot: Style 2. With mean Lhat(t) and absolute error

% A figure with four curves  
sim_index = 1;
figure(2)
plot(vec_time,mat_llr_exact(:,sim_index),'b-', ...
     vec_time,mat_llr_appro(:,sim_index),'r-', ...
     tv_ia,yv_ia(:,2),'m-', ... 
     vec_time,mean_abs_error,'k', ...
     'Linewidth',2);  
legend({'L(t)','$\hat{L}(t)$','mean $\hat{L}(t)$','mean $|L(t)-\hat{L}(t)|$'}, ...
        'Location','East', ...
        'FontWeight','Bold','FontSize',20,'Interpreter','latex')
xlabel('time','FontSize',20) 
ylabel('Log-likelihood Ratio','FontSize',20)
eval(['print -depsc plot_ia_',int2str(test_num)])
 
