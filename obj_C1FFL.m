function y_stacked = objC1FFL_AD(x0,time_vec,amp_set,dur_set,para_x,input_base)
%
% This function solves the differential equations (DEs) for C1-FFL 
% for inputs which are rectangular pulses. Each rectangular pulse
% is defined by an amplitude and a duration. 
% 
% The input amp_set specifies a set of amplitudes
% The input dur_set specifies a set of durations
% 
% This functions solves the DE by using inputs coming from the 
% Cartesian product of amp_set and dur_set 
%
% Inputs:
%   x0              Contains the parameters of the Y and Z part for C1-FFL 
%   time_vec        Time vector for ODE45 
%   amp_set         Set of input amplitudes
%   dur_set         Set of input duration 
%   para_x          Parameter of the X part of C1-FFL. These parameters 
%                   are not used for optimisation. 
%   input_base      A vector with 3 elements to be used with 
%                   ode_C1FFL 
%                   [off_amplitude on_amplitude duration] 
%                   This is basically used to pass the value of
%                   off_amplitude to this function
%
% Outputs:
%   y_stacked       The outputs of the C1-FFL stacked together 

% Extract the parameters of the Z and Y part from x0
para_z = [x0(1:3) 0];
para_y = [1 x0(4:7)]; 

% To store the results in the simulation of C1-FFL in a 
% 3-dimensional array
y_amp_dur = zeros(length(time_vec),length(amp_set),length(dur_set)); 

% Initial condition for DEs for C1-FFL
init_ode = zeros(3,1); 

% Loop through all the amplitudes and durations 
for aa = 1:length(amp_set)
    for dd = 1:length(dur_set)
        input_para = input_base; % [off_amplitude on_amplitude duration]
        input_para(2) = amp_set(aa);
        input_para(3) = dur_set(dd); 
        
        [dummy,yv_tmp] = ode45(@(t,x) ode_C1FFL(t,x,input_para,para_x, para_y, para_z),time_vec,init_ode);
        y_amp_dur(:,aa,dd) = yv_tmp(:,3); 
    end
end

% Turn the 3-D array into a vector 
y_stacked = y_amp_dur(:); 