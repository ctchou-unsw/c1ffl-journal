function dx = ode_ia(t, x, para_input, para_sys, para_ref, input_type) 
% 
% This function defines the system of ODEs that models 
% (1) Activation and deactivation of x
% (2) Approximate computation of log-likelihood ratio, i.e. the 
%     intermediate approximation in the paper 
% 
% The function is to be used with an ODE solver, e.g. ode45 
%
% Inputs:
%       t               time 
%       x               state of the system (2 elements) 
%                       x(1) = x_star(t)  
%                       x(2) = \hat{L}(t) 
%       para_input      (vector) [off_amplitude on_amplitude duration]
%       para_sys        (vector) [kx dx Mx]
%       para_ref        Parameters of the reference system
%
% Output:
%       dx              dx/dt at time t 
% 
% Chun Tung Chou, UNSW 6/1/18, 16/8/18
%  

    % Default input type
    if nargin < 6
        input_type = 'rect';
    end
    
    % Unpack para_sys which contains the system parameters 
    % para_sys = [kx dx Mx];
    k_x = para_sys(1);
    d_x = para_sys(2);
    M_x = para_sys(3);
      
    % Unpack para_ref 
    % input_ref = [basel_conc sig_conc dur_short dur_long];    
    amp_ref_0 = para_ref(1);
    amp_ref_1 = para_ref(2);   
    dur_short_ref = para_ref(3);
    dur_long_ref = para_ref(4);       
            
    % Initialise dx
    dx = zeros(size(x)); 
   
    % Calculate the input  
    switch input_type
        case 'rect'
            % Unpack para_input 
            % The format is: 
            % (vector) [off_amplitude on_amplitude duration]
            amp0 = para_input(1);
            amp1 = para_input(2);
            dur = para_input(3);
            if (t >= 0) && (t <= dur)        
                u = amp1;
            else
                u = amp0;
            end  
        case 'tri'
            % Unpack para_input 
            % The format is: 
            % [d1_tri d2_tri amp_tru] 
            % Ramp up from 0 to amp_tru in [0,d1_tri]
            % Ramp down from amp_tri to 0 in [d1_tri,d2_tri]
            d1_tri = para_input(1);
            d2_tri = para_input(2);
            amp_tri = para_input(3);
            % Calculate the input signal 
            if (t >= 0) && (t <= d1_tri)
                u = amp_tri*(t/d1_tri);
            elseif (t > d1_tri) && (t <= d2_tri)
                u = amp_tri*(1-(t-d1_tri)/(d2_tri-d1_tri));
            else
                u = 0; 
            end    
    end 
    
    % Calculate the weighting function 
    if (t >= 0) && (t <= dur_short_ref)        
        r0 = amp_ref_1;
    else
        r0 = amp_ref_0;
    end     

    if (t >= 0) && (t <= dur_long_ref)        
        r1 = amp_ref_1;
    else
        r1 = amp_ref_0;
    end      

    log_term = log(r1/r0);
    diff_term = (r1-r0);

    % The derivatives 
    dx(1) = k_x * u * (M_x - x(1)) - d_x * x(1); % Amount of X* 
    % 
    % The matching term
    matching = max(log_term - diff_term/u,0);
    dx(2) = d_x * x(1) * matching; 

end

