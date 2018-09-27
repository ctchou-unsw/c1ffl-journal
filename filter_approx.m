function llr_approx = filter_approx(time_vec,xstar_vec,input_para,input_ref,para_sys)
%
%  Inputs:
%   time_vec    time vector from SSA simulation 
%   xstart_vec  the number of X* molecules (from SSA) 
%   input_para  (vector) [off_amplitude on_amplitude duration]
%   input_ref   (vector) [off_amplitude_ref on_amplitude_ref duration_short_ref duration_long_ref]
%   para_sys    (vector) [kx dx Mx]
% 
%  Outputs:
%   llr_approx  The approximate LLR 
%
%  Chun Tung Chou, UNSW. This version 16/8/18  
% 

% The rate constant for deactivation reaction
dx = para_sys(2);

% The parameters of the input signal 
% (vector) [off_amplitude on_amplitude duration]
dur_input = input_para(3);
on_amp_input = input_para(2);
off_amp_input = input_para(1);

% The parameters of the reference signal 
% (vector) [off_amplitude_ref on_amplitude_ref duration_short_ref duration_long_ref]
dur_long_fef = input_ref(4);
dur_short_ref = input_ref(3); 
on_amp_ref = input_ref(2);

% Terms from the right hand side of the ODE 
% pi(t) 
pi_t = ones(size(time_vec));
pi_t(time_vec >= dur_long_fef) = 0;
pi_t(time_vec <  dur_short_ref) = 0;

% The input signal and the weights 
u_t = zeros(size(time_vec));  
u_t(time_vec < dur_input) = on_amp_input;
u_t(time_vec >= dur_input) = off_amp_input;

weights_t = max(log(on_amp_ref/off_amp_input) - (on_amp_ref-off_amp_input) ./ u_t,0);

% RHS of the ODE 
rhs = dx * xstar_vec .* pi_t .* weights_t;  

llr_approx = cumtrapz(time_vec,rhs); 

