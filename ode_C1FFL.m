function dx = ode_C1FFL(t, x, para_input, para_x, para_y, para_z, input_type) 
% 
% This function is to be used with ode45 to solve the differential
% equations (DEs) describing an C1-FFL 
%
% The DEs are 
%   dx/dt = -k_x (M_x - x) + d_x x 
%   dy/dt = k_y * Hill_y(x) - d_y y 
%   dz/dt = x Hill_z(y) - d_z z 
%
% where 
%   Hill_y(x) = h_y * x^n_y / (x^n_y + Kd_y^n_y)
%   Hill_z(y) = h_z * y^n_z / (y^n_z + Kd_z^n_z)
%
% Inputs:
%       para_input      Parameters of the input which is assumed to be 
%                       a rectangular pulse. The pulse parameters are 
%                       [off_amplitude on_amplitude duration]
%       para_x          Parameters of the DE for x
%                       [k_x d_x M_x]
%       para_y          Parameters of Hill_y(x) in the DE for dy/dt 
%                       [k_y h_y Kd_y n_y d_y]
%       para_z          Parameters of Hill_z(y) in the DE for dz/dt 
%                       [h_z n_z Kd_z d_z]                 
% 
% Chun Tung Chou, UNSW. 18/8/18

    % Default input type
    if nargin < 7
        input_type = 'rect';
    end
    
    % Unpack para_x which is the vector  
    % [k_x d_x M_x]
    k_x = para_x(1);
    d_x = para_x(2);
    M_x = para_x(3);

    % Unpack para_y which is the vector [k_y h_y Kd_y n_y d_y]
    k_y = para_y(1);
    h_y = para_y(2);
    Kd_y = para_y(3);
    n_y = para_y(4);
    d_y = para_y(5);    
        
    % Unpack para_z which is the vector [h_z n_z Kd_z d_z]
    h_z = para_z(1);     
    n_z = para_z(2);
    Kd_z = para_z(3);
    d_z = para_z(4); 
    
    % Initialise dx
    dx = zeros(3,1); 
   
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
    
    % The derivatives 
    dx(1) = k_x * u * (M_x - x(1)) - d_x * x(1); % Amount of X* 
    dx(2) = k_y * h_y * x(1)^n_y / (x(1)^n_y + Kd_y^n_y)- d_y * x(2);
    dx(3) = h_z * x(1) * x(2)^n_z / (x(2)^n_z + Kd_z^n_z) - d_z * x(3);     
end

