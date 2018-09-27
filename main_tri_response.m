% This script compares the response of the LLR detector with that of a 
% C1-FFL when the input is a triangular pulse 

% Load the optimised parameters 
load c1ffl_optpara 

% Simulation end time 
time_end = 2*dur_long;
time_span = [0 time_end];

% Specifying the inputs of the triangular pulse 
inputTri = [40 80 4*on_amp]; % sigConc 10.71

% Solve the ODE for the detection system  
init_ia = zeros(2,1);

[tv_ia,yv_ia] = ode45(@(t,x) ode_ia(t,x,inputTri,para_sys,input_ref,'tri'),time_span,init_ia);

% Simulate the C1FFL 
init_ffl = zeros(3,1);
[tv_ODE_FFm4z,yv_ODE_FFm4z] = ode45(@(t,x) ode_C1FFL(t,x,inputTri ,para_sys, para_y_opt, para_z_opt,'tri'),time_span,init_ffl);

figure(2)
plot(tv_ia,yv_ia(:,2),'b',tv_ODE_FFm4z,yv_ODE_FFm4z(:,3),'r','linewidth',3)
legend({'$\hat{L}(t)$','C1-FFL'},'Location','NorthWest','Interpreter','latex','Fontsize',14)
ylabel('$\hat{L}(t)$ / C1-FFL output','Interpreter','latex','Fontsize',14)
xlabel('time','Interpreter','latex','Fontsize',14)
print -depsc plot_tri_resp
% title('Triangular input s(t)')



