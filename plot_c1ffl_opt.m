% Load the optimised parameters 
load c1ffl_optpara % 

% Time vector for computation 
time_vec = 0:0.5:100;
% The different durations 
dur_range_2 = [20 50 90];
% Time for amplitdue 
time_check = 70;

amp_at_end_of_pulse_ia_check = zeros(length(amp_set),1);
amp_at_end_of_pulse_c1ffl_check = zeros(length(amp_set),1);

mat_y_llr = zeros(length(time_vec),length(amp_set),length(dur_range_2)); 
mat_y_opt = zeros(length(time_vec),length(amp_set),length(dur_range_2)); 

% Simulate the intermediate approximation   
for aa = 1:length(amp_set)
    for dd = 1:length(dur_range_2)
        input_para = input_base; % [basalConc sigConc durLong];
        input_para(2) = amp_set(aa);
        input_para(3) = dur_range_2(dd); 
        
        [dummy,yv_tmp] = ode45(@(t,x) ode_ia(t,x,input_para,para_sys,input_ref),time_vec,init_ia);
        mat_y_llr(:,aa,dd) = yv_tmp(:,2); 
        amp_at_end_of_pulse_ia_check(aa) = yv_tmp(time_check,2);
    end
end

% The optimised system 
for aa = 1:length(amp_set)
    for dd = 1:length(dur_range_2)
        input_para = input_base; 
        input_para(2) = amp_set(aa);
        input_para(3) = dur_range_2(dd); 
        
        % Simulate optimised system 
        [dummy,yv_tmp] = ode45(@(t,x) ode_C1FFL(t,x,input_para, para_sys, para_y_opt, para_z_opt),time_vec,init_ffl);
        mat_y_opt(:,aa,dd) = yv_tmp(:,3); 
        amp_at_end_of_pulse_c1ffl_check(aa) = yv_tmp(time_check,3);                      
    end
end

%% Result 1: A chosen amplitude but different durations 
aa = 4; % amplitude = 10.7143
figure(1)
h = plot(time_vec,mat_y_llr(:,aa,3),'b: ',time_vec,mat_y_opt(:,aa,3),'r: ', ...
         time_vec,mat_y_llr(:,aa,2),'b--',time_vec,mat_y_opt(:,aa,2),'r--', ...
         time_vec,mat_y_llr(:,aa,1),'b-' ,time_vec,mat_y_opt(:,aa,1),'r-',  ...
         'Linewidth',3); 
legend({['$\hat{L}(t)$, d = ',int2str(dur_range_2(3))],['C1-FFL, d = ',int2str(dur_range_2(3))], ...
        ['$\hat{L}(t)$, d = ',int2str(dur_range_2(2))],['C1-FFL, d = ',int2str(dur_range_2(2))], ...
        ['$\hat{L}(t)$, d = ',int2str(dur_range_2(1))],['C1-FFL, d = ',int2str(dur_range_2(1))]}, ...
       'Location','NorthWest','FontWeight','Bold','FontSize',20,'Interpreter','latex')
xlabel('time','FontSize',20) 
ylabel('$\hat{L}(t)$ / C1FFL output','FontSize',20,'Interpreter','latex')
xticks([0 40 80])
yticks(0:500:1000) 
ax = ancestor(h, 'axes');
yrule = ax{1}.YAxis;
yrule.FontSize = 18;
xrule = ax{1}.XAxis;
xrule.FontSize = 18;
print -depsc plot_Lhat_C1FFL_fixedA

%% Result 1_b: A chosen amplitude but different durations 
aa = 15; % amplitude = 40.1786
figure(2)
h = plot(time_vec,mat_y_llr(:,aa,3),'b: ',time_vec,mat_y_opt(:,aa,3),'r: ', ...
         time_vec,mat_y_llr(:,aa,2),'b--',time_vec,mat_y_opt(:,aa,2),'r--', ...
         time_vec,mat_y_llr(:,aa,1),'b-' ,time_vec,mat_y_opt(:,aa,1),'r-',  ...
         'Linewidth',3); 
legend({['$\hat{L}(t)$, d = ',int2str(dur_range_2(3))],['C1-FFL, d = ',int2str(dur_range_2(3))], ...
        ['$\hat{L}(t)$, d = ',int2str(dur_range_2(2))],['C1-FFL, d = ',int2str(dur_range_2(2))], ...
        ['$\hat{L}(t)$, d = ',int2str(dur_range_2(1))],['C1-FFL, d = ',int2str(dur_range_2(1))]}, ...
       'Location','NorthWest','FontWeight','Bold','FontSize',20,'Interpreter','latex')
xlabel('time','FontSize',20) 
ylabel('$\hat{L}(t)$ / C1FFL output','FontSize',20,'Interpreter','latex')
xticks([0 50 100])
yticks(0:5000:10000) 
ax = ancestor(h, 'axes');
yrule = ax{1}.YAxis;
yrule.FontSize = 18;
xrule = ax{1}.XAxis;
xrule.FontSize = 18;
print -depsc plot_Lhat_C1FFL_fixedA_2


%% Result 2: A chosen time point but different input amplitudes

figure(3)
h = ...
plot(amp_set,amp_at_end_of_pulse_ia_check,'b-', ...
     amp_set,amp_at_end_of_pulse_c1ffl_check ,'r-','Linewidth',3);
legend({'$\hat{L}$','C1-FFL'}, ...
        'Location','NorthWest','FontWeight','Bold','FontSize',20,'Interpreter','latex')
xlabel('Input amplitude $a$','FontSize',20,'Interpreter','latex') 
ylabel('$\hat{L}$ / C1FFL output','FontSize',20,'Interpreter','latex')
xticks(0:20:40)
yticks(0:2000:7000) 
ax = ancestor(h, 'axes');
yrule = ax{1}.YAxis;
yrule.FontSize = 18;
xrule = ax{1}.XAxis;
xrule.FontSize = 18;
print -depsc plot_Lhat_C1FFL_fixedTime