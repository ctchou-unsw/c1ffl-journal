% This script generates the data for Figure 3 for illustrating the 
% behaviour of the intermediate approximation 
% 
% The aims are to compare: 
% (1) Calculation of L(t) over one realisation
% (2) Calculation of \hat{L}(t) over one realisation
% (3) Comparing the error |L(t) - \hat{L}(t)| over 100 realisations
% (4) Calculation of the mean of \hat{L}(t)
% 
% Chun Tung Chou, UNSW 
% 

%% Problem parameters 
% Parameters for the X-X* reaction cycle 
kx = 0.02;  % Reaction rate constants S + X -> S + X* 
dx = 0.5;   % Reaction rate constants X* -> X
Mx = 100;   % Total of X and X* 

% Input is defined by 3 parameters: 
% OFF (basal) amplitude, ON amplitude, duration of the signal
% [basel amplitude duration]
% 
% The concentrations are determined by ratio of the steady state number of
% X* to the total number of X*
ratio_basal = 0.01; % steady state X* to Mx for basal concentration 
ratio_sig = 0.3;    % steady state X* to Mx for signal concentration 
dur_short = 5; 
% dur_long = 60; % Select at test configuration 

% Duration of the reference signal 
dur_short_ref = 10;
dur_long_ref = 60; 

% Deduce the concentrations to reach the steady state concentration
% specified by ratioBasal and ratioSig 
basal_conc = (dx/kx)/(1/ratio_basal-1);  % 0.25 
sig_ref_conc = (dx/kx)/(1/ratio_sig-1);  % 10.7 

% The parameters sig_conc and dur_long are specified in the test
% configuration in the next cell

%% Choose the test configuration
test_num = 3;
switch test_num
    case 1  % Sub-figure (a) and (b) 
        sig_conc = sig_ref_conc;   % ON-amplitdue for input 
        dur_long = 70;             
    case 2  % Sub-figure (c)  
        sig_conc = sig_ref_conc;
        dur_long = 50;
    case 3  % Sub-figure (d)
        sig_conc = 37.5;
        dur_long = 40; 
end     

%% Parameter vectors 
input_short = [basal_conc sig_conc dur_short];
input_long = [basal_conc sig_conc dur_long];
input_ref = [basal_conc sig_ref_conc dur_short_ref dur_long_ref];
para_sys = [kx dx Mx]; 

%% Simulation end time 
time_end = 2*dur_long;
time_span = [0 time_end];

%%  Intermediate approximation - mean of \hat{L}(t) 
% Parameters for the intermediate approximation ODE 
% ODE of reaction-cycle and intermediate approximation (IA) 
init_ia = zeros(2,1);
[tv_ia,yv_ia] = ode45(@(t,x) ode_ia(t,x,input_long,para_sys,input_ref),time_span,init_ia);

%% SSA simulations
% Two sets of simulation 
% (1) log-likelihood calculation 
% (2) intermediate approximation 

% Simulation parameters 
n_sim = 100; % Simulate SSA 100 times
vec_time = 0:0.1:time_end;

% Exact Log-likelihood calculation 
tv_ref = 0:0.1:time_end;
amp_ref_short = repmat(basal_conc,1,length(tv_ref));
amp_ref_short(tv_ref <= dur_short_ref) = sig_ref_conc;
amp_ref_long = repmat(basal_conc,1,length(tv_ref));
amp_ref_long(tv_ref <= dur_long_ref) = sig_ref_conc;

% Storage the result for log-likelihood computation 
mat_llr_exact = zeros(length(vec_time),n_sim); % Exact  
mat_llr_appro = zeros(length(vec_time),n_sim); % Interemediate approximation 

% Simulate and calculate 
for i = 1:n_sim
    % Simulate the reaction cycle 
    [tv_ssa_x_long,yv_ssa_x_long] = ssa_simple_cycle(para_sys,input_long,time_end);
    
    % Exact likelihood ratio computation     
    [tv_llr_exact,llr_exact] = filter_exact(tv_ssa_x_long,yv_ssa_x_long(:,2),tv_ref,input_ref,para_sys);
    llr_exact_uniform_time = interp1(tv_llr_exact,llr_exact,vec_time);
    mat_llr_exact(:,i) = llr_exact_uniform_time;
    
    % Intermediate approximation 
    llr_appro = filter_approx(tv_ssa_x_long,yv_ssa_x_long(:,2),input_long,input_ref,para_sys);
    llr_appro_uniform_time = interp1([tv_ssa_x_long ; time_end],[llr_appro ; llr_appro(end)],vec_time);
    mat_llr_appro(:,i) = llr_appro_uniform_time; 
end

% Compute the mean absolute error
mean_abs_error = mean(abs(mat_llr_exact-mat_llr_appro),2);

%% Plots

% Figure 1: One simulation from each 
sim_index = 1;
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

eval(['save eval_int_approx_v1_data_ex',int2str(test_num)])