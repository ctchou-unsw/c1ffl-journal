function [t_ssa1,y_ssa1] = ssa_simple_cycle(vec_para_sys,vec_input,time_end)
%
% Performs SSA simulation of an activation-deactivation cycle  
%
%      S + X  -> S + X* at a rate of kx
%          X* -> X      at a rate of dx
%
% The total number of X and X* is a constant Mx 
% 
% Inputs:
%   vec_para_sys    (vector) [kx dx Mx]
%   vec_input       (vector) [off_amplitude on_amplitude duration]
%   time_end        (scalar) end time of SSA simulation 
%
% Outputs:
%   t_ssa1:  vector of reaction time instances 
%   y_ssa1:  vector with 2 columns. First column = #X. Second column = #X* 
%
% Chun Tung Chou, UNSW
% 

%% Unpack the system parameters from the vector vec_para_sys 
% Format: vec_para_sys = [kx dx Mx];
kx = vec_para_sys(1);
dx = vec_para_sys(2);
Mx  = vec_para_sys(3);


%% Unpack the input parameter from vec_input 
% The format is: 
% [off_amplitude on_amplitude duration]
amp0 = vec_input(1);
amp1 = vec_input(2);
dur = vec_input(3);

%% Prepare two vectors for SSA 
% tDetInput = time instances at which the input changes
% inputPromotor = the input function for the promotor
time_input = [0 dur];
amp_input = [amp1 amp0];
% Find the number of time segments  
num_segments = length(time_input);    

% Create the stoichiometric matrix and structure needed for SSA
num_reactions = 2;
re_activation = 1;
re_deactivation = 2;
% 
num_species = 2; % Unbound Promotor, Bound Promotor, 
chem_unbound = 1;
chem_bound = 2;
% The stochiometric matrix 
mat_stoi = zeros(num_reactions,num_species);
mat_stoi(re_activation   ,[chem_unbound chem_bound]) = [-1  1]; % Bound reaction
mat_stoi(re_deactivation ,[chem_unbound chem_bound]) = [ 1 -1]; % Unbound reaction
% The 'rcv' vector - for the reaction rate constants
vec_rrc = zeros(1,num_reactions);
vec_rrc(re_activation   ) = kx*amp1;  % This will be overwritten later 
vec_rrc(re_deactivation ) = dx;
% The 'rsi' vector - for the state in the reactions
vec_rsi = zeros(1,num_reactions);
vec_rsi(re_activation   ) = chem_unbound;
vec_rsi(re_deactivation ) = chem_bound;
ssaPara = struct('rcv',vec_rrc,'rsi',vec_rsi);

% Initialisation for SSA 
t_ssa0 = 0;         % Simulation time interval is [t_ssa0,t_ssa1]
% Initial vector 
x0 = zeros(size(mat_stoi,2),1); 
x0(chem_unbound) = Mx;   

% These are the outputs 
y_ssa = [];
t_ssa = [];

% function handle for propensity function 
pfun = @propensity_fun;

% Perform SSA for each interval where there are no input emissions in
% between 
for i = 1:num_segments      % Loop through all non-zero input 
    % This is to determine the upper limit of simulation interval
    % If it is inbetween two input, use the next input time as the end of
    % simulation interval, otherwise use the end time 
    if i < num_segments
        t_ssa1 = time_input(i+1); 
    else % Last interval 
        t_ssa1 = time_end;
    end
    
    % Adjust the reaction rate constants according to the input   
    ssaPara.rcv(re_activation) = kx * amp_input(i); 
    
    % Specify time span 
    tspan = [t_ssa0 t_ssa1]; 
    % Run SSA 
    [t_ssa_tmp,x_ssa_tmp] = firstReactionMethod(mat_stoi, pfun, tspan, x0, ssaPara); 
    % x_ssa_tmp has the same number of columns as the number of species     
    % 
    % Concatenate the outputs 
    t_ssa = [t_ssa ; t_ssa_tmp];
    y_ssa = [y_ssa ; x_ssa_tmp]; % # signalling molecules in receiver voxel
    % update for the next non-zero input - except for the last round
    % Add molecules to the transmitter voxel 
    if i < num_segments 
        t_ssa0 = t_ssa1;
        x0 = x_ssa_tmp(end,:);
        x0 = x0(:);
    end
end

% Remove the duplications due to segment changes 
index_y_change = [1 ; 1+find(sum(abs(diff(y_ssa)),2) ~= 0)];
t_ssa1 = t_ssa(index_y_change);
y_ssa1 = y_ssa(index_y_change,:);


end % end of function 


% compute the reaction rate 
function rrate = propensity_fun(x,para)
%
% x is a row vector 
% rrate is a column vector 

if any(x < 0) 
    error('negative state')
end 

% unpack parameters 
rcv = para.rcv; % reaction rate 
rsi = para.rsi; % state relevant to a state 

rrate = rcv .* x(rsi);
rrate = rrate(:); 

% This is for simulate recetpor binding whose rate is a product of two states 
% Note: This method requires co-ordination with a correct initial condition
% 
if isfield(para,'two')  % if a two species reaction is defined 
    % #rows specifies the number of such reactions 
    % each row has 4 elements 
    % [reaction_number state1 state2 kinetic_constant]
    nRTR = size(para.two,1); 
    for i = 1:nRTR
        % Extract i-th row to cRe (current Reaction)
        cRe = para.two(i,:);
        reactionIndex = cRe(1);
        kineticConstant = cRe(4);
        numSpecies1 = x(cRe(2));
        numSpecies2 = x(cRe(3));
        rrate(reactionIndex) = kineticConstant * numSpecies1 * ...
                               numSpecies2;
    end 
end    




end