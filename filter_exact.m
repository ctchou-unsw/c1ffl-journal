function [tv_all,llr] = filter_exact(tv_ssa,nxstar_ssa,tv_uni,para_ref,para_sys)
% 
% To compute log-likelihood ratio (LLR) from the SSA simulation data 
% 
% Inputs
%   tv_ssa      (vector) time vector from SSA
%   nxstar_ssa  (vector) number of X* molecules 
%   tv_uni      (vector) A uniform time vector for computing LLR 
%   para_ref    (vector) Parameters of the reference signals  
%               [off_amplitude_ref on_amplitude_ref duration_short_ref duration_long_ref]
%   para_sys    (vector) [kx dx Mx]
% 
%
% Outputs:
%   tv_all      (vector) time vector for LLR 
%   llr         (vector) vector of LLR values
%
%
% Methods: Trapezoidal rule for the continuous part of integration using
%          cumtrapz()

% Extract parameters from para_sys 
kx = para_sys(1);
Mx = para_sys(3);

% The reference amplitudes and durations
amp_ref0 = para_ref(1);
amp_ref1 = para_ref(2);
dur_short_ref = para_ref(3);
dur_long_ref = para_ref(4);   

% Extrapolate tv_ssa and nxstar_ssa to include the end time 
tv_uni = tv_uni(:); 
tv_ssa = [tv_ssa ; tv_uni(end)];
nxstar_ssa = [nxstar_ssa ; nxstar_ssa(end)];

% Merge the time vectors tvc and tvw 
tv_all = sort([tv_ssa ; tv_uni(2:end-1)]); 

% Determine the times at which activations occur 
nxstar_all = interp1(tv_ssa,nxstar_ssa,tv_all,'previous'); 
nxstar_diff = [0 ; diff(nxstar_all)];
nxstar_activate = nxstar_diff; 
nxstar_activate(nxstar_activate == -1) = 0; 

% Determine the function pi(t) 
pi_t = zeros(size(tv_all));
pi_t((tv_all > dur_short_ref) & (tv_all <= dur_long_ref)) = 1; 

% Use cumsum and cumtrapz to do the integration 
llr_term1 = cumsum(pi_t .* nxstar_activate) * log(amp_ref1/amp_ref0);
llr_term2 = kx*(amp_ref1-amp_ref0)*cumtrapz(tv_all,pi_t .* (Mx-nxstar_all));

% Substraction to get LLR 
llr = llr_term1 - llr_term2; 