% This scripts determines the parameters of the C1FFL   
% that fit a particular intermediate approximation 

% The activation-deactivation cycle parameters
kx = 0.02;  % Reaction rate constants for X*
dx = 0.5; 
Mx = 100;   % Total amount of X and X* 

% Input is defined by 3 parameters: 
% OFF-amplitude (= basal), ON-amplitude, signal duration 
% [OFF-amplitude, ON-amplitude, duration]
% The concentrations are determined by ratio of the steady state number of
% X* to the total number of X*
ratio_basal = 0.01; % steady state X* to Mx for basal concentration 
ratio_on = 0.3;     % steady state X* to Mx for signal concentration 
dur_short = 10; 
dur_long = 80; 
dur_short_ref = 10;
dur_long_ref = 80; 

% Simulation end time 
time_end = 2*dur_long;
time_span = [0 time_end];

% Deduce the concentrations to reach the steady state concentration
% specified by ratioBasal and ratioSig 
off_amp = (dx/kx)/(1/ratio_basal-1);
on_amp  = (dx/kx)/(1/ratio_on-1);

%%  Assemble parameters into vectors 
input_short = [off_amp on_amp dur_short];
input_long = [off_amp on_amp dur_long];
input_ref = [off_amp on_amp dur_short_ref dur_long_ref];
para_sys = [kx dx Mx]; 


%% Fit the matching function to a Hill function
% The matching function is:
% log(on_amp/off_amp) - (on_amp - off_amp) / x 
% 
% This is so that we can use the fitted Hill function 
% to come out with the initial condition for one of the
% Hill functions 
% 
match1 = log(on_amp/off_amp);
match2 = on_amp - off_amp;
% Create the independent and dependent vectors 
step_size = 0.1;
lower_fitting_limit = ceil(match2/match1/step_size)*step_size;
upper_fitting_limit = 20;
vec_ind = lower_fitting_limit:step_size:upper_fitting_limit; 
vec_dep = match1 - match2 ./ vec_ind;
curve_to_be_fitted = fittype('k*x^n/(Kd^n+x^n)','independent','x');
coeffnames(curve_to_be_fitted); % Kd, k, n 
n_guess = 1.5; Kd_guess = 2; h_guess = 3.5; % Initial guess 
vec_dep_guess = h_guess*vec_ind.^n_guess./(vec_ind.^n_guess+Kd_guess^n_guess); 
fit_output = fit(vec_ind(:),vec_dep(:),curve_to_be_fitted, ...
                'StartPoint',[Kd_guess h_guess n_guess]);
vec_dep_fit = fit_output(vec_ind(:));
h_m = fit_output.k;
Kd_m = fit_output.Kd;
n_m = fit_output.n;

figure(10)
plot(vec_ind,vec_dep,vec_ind,vec_dep_guess,vec_ind,vec_dep_fit)
legend('matching','Init','Opt')
title('Fit of the matching function') 

%% Fit to the coherent Type-1 FFL 
%
% The aim of this part is to fit a C1-FFL to the intermediate 
% approximation 


% The DEs are C1-FFL are: 
%   dx/dt = -k_x (M_x - x) + d_x x 
%   dy/dt = k_y * Hill_y(x) - d_y y 
%   dz/dt = x Hill_z(y) - d_z z 
%
% where 
%   Hill_y(x) = h_y * x^n_y / (x^n_y + Kd_y^n_y)
%   Hill_z(y) = h_z * y^n_z / (y^n_z + Kd_z^n_z)
% 
% Initial parameters for the C1-FFL 
% para_y0         Parameters of Hill_y(x) in the DE for dy/dt 
%                 [k_y h_y Kd_y n_y d_y]
% para_z          Parameters of Hill_z(y) in the DE for dz/dt 
%                 [h_z n_z Kd_z d_z]  
%  
para_y0 = [0.05 h_m Kd_m n_m 0.10];  
para_z0 = [2   5    1.3 0];       
x0_init = [para_z0(1:3) prod(para_y0([1 2])) para_y0(3:5)];  

% initial condition for the ODEs  
init_ffl = zeros(3,1); 
init_ia = zeros(2,1);

%% Optimisation of the parameters in the Y and Z parts of c1-FFL
%
% In the optimisation, we will use rectangular pulse as an input to
% the intermediate approximation. We do that for rectangular pulses whose
% amplitude and duration come from, respectively, the sets amp_set and
% dur_set. This is the data to be fittted to.
%
% We determine the C1-FFL parameters to fit this data. 
% 

% Define the amplitude and duration set 
amp_set = (0.25:0.25:4)*on_amp;
dur_set = 70;  % Assuming to be a single value at the moment
time_end_opt = dur_set+10; 
time_vec = 0:0.5:time_end_opt; 
input_base = [off_amp 0 0]; % Needed for obj_C1FFL

% Storage the time profiles
% The data for fitting
y_fitting_array = zeros(length(time_vec),length(amp_set),length(dur_set)); 
% The array for optimisation 
y_opt_array = zeros(length(time_vec),length(amp_set),length(dur_set)); 
% The C1FFL output for the initial condition in optimisation
y_init_array = zeros(length(time_vec),length(amp_set),length(dur_set));
% Time at the end of the pulse (note: assuming singleton dur_set)
time_end_of_pulse = find(time_vec == dur_set);

y_at_pulse_end_ia = zeros(length(amp_set),1);
y_at_pulse_end_opt = zeros(length(amp_set),1);
y_at_pulse_end_init = zeros(length(amp_set),1);

for aa = 1:length(amp_set)
    for dd = 1:length(dur_set)
        input_para = input_base; % [off_amp on_amp durLong];
        input_para(2) = amp_set(aa);
        input_para(3) = dur_set(dd); 
        
        [dummy,yv_tmp] = ode45(@(t,x) ode_ia(t,x,input_para,para_sys,input_ref),time_vec,init_ia);
        y_fitting_array(:,aa,dd) = yv_tmp(:,2); 
        y_at_pulse_end_ia(aa) = yv_tmp(time_end_of_pulse,2);
    end
end
ydata = y_fitting_array(:); 


% Optimisation options 
optim_options = optimset('Display','iter','MaxFunEvals',1000); %1000
% Lower bounds
lb = zeros(1,5); lb(4) = 1.01;
% Optimisation using lsqcurvefit
x1_opt = lsqcurvefit(@(x0,xdata) obj_C1FFL(x0,xdata,amp_set,dur_set,para_sys,input_base), ...
                     x0_init,time_vec,ydata,lb,[],optim_options); 
 
% Plot the amplitudes 
para_y_opt = [1 x1_opt(4:7)];   
para_z_opt = [x1_opt(1:3) 0]; 

% Simulate the C1-FFL with optimised parameters 
for aa = 1:length(amp_set)
    for dd = 1:length(dur_set)
        input_para = input_base; % [off_amp on_amp durLong];
        input_para(2) = amp_set(aa);
        input_para(3) = dur_set(dd); 
        
        % Simulate optimised system 
        [dummy,yv_tmp] = ode45(@(t,x) ode_C1FFL(t,x,input_para ,para_sys, para_y_opt, para_z_opt),time_vec,init_ffl);
        y_opt_array(:,aa,dd) = yv_tmp(:,3); 
        y_at_pulse_end_opt(aa) = yv_tmp(time_end_of_pulse,3);
               
        % Simulate initial system 
        [dummy,yv_tmp] = ode45(@(t,x) ode_C1FFL(t,x,input_para ,para_sys, para_y0, para_z0),time_vec,init_ffl);
        y_init_array(:,aa,dd) = yv_tmp(:,3);
        y_at_pulse_end_init(aa) = yv_tmp(time_end_of_pulse,3);

        
    end
end

%% Plot
% For each amplitude 
for aa = 1:length(amp_set)
    figure(10+aa)
    plot(time_vec,y_fitting_array(:,aa,1),'r', ...
         time_vec,y_opt_array(:,aa,1),'b',time_vec,y_init_array(:,aa,1),'m')
    legend('True','Opt','Init','Location','SouthEast')
    title(['Amp = ',num2str(amp_set(aa))])    
end 

% For the end of pulse 
figure
plot(amp_set,y_at_pulse_end_ia,'--',amp_set,y_at_pulse_end_opt,'-')
xlabel('Pulse ON-amplitude')
ylabel('Output at end of the pulse')

% Earlier 
%      0          8     4.63273e+09                      1.37e+10
%      1         16     4.63273e+09        2.89715       1.37e+10      
%      2         24     5.65375e+08       0.724287       2.45e+09      
%      3         32     3.16022e+08        1.44857       1.49e+08      
%      4         40     3.16022e+08        2.89715       1.49e+08   
%      5         48     2.86447e+08       0.724287       2.54e+08  

% Now
% x0_init
% 2.0000    5.0000    1.3000    0.1611    5.1714    2.7240    0.1000
% Iteration 
%      0          8     4.62448e+09                      1.37e+10