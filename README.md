%% file_info.m 
% 
% This file summarises the information related to Matlab files 
% that were used to generate the results in the paper:
%
% Title: Detection of persistent signals and its relation to coherent 
% feedforward loops
% Author: Chun Tung Chou, University of New South Wales, Sydney, Australia
% E-mail: c.t.chou@unsw.edu.au
%
% This version: 26/09/2018 
 
%% File for Figure 5
% For the intermediate approximation

% Main file: main_int_approx

% Auxiliary files: 
%   ode_ia              Use with ODE45 to solve the intermediate approximation 
%   ssa_simple_cycle    SSA simulation of the reaction cycle 
%   filter_exat         The filter for computint the exact LLR 
%   filter_approx       The filter for computing the intermediate approximation of LLR
%   plot_int_approx     For plotting the graphs 
%
%   Note that the file ssa_simple_cycle uses the m-file
%   firstReactionMethod, which implements Gillespie's SSA algorithm.
%   The file is written by Nezar Abdennur and can be obtained from
%   https://au.mathworks.com/matlabcentral/fileexchange/34707-gillespie-stochastic-simulation-algorithm

%% File for Figure 6 

% Plot the results in Figure 6a, b, c
% Main file: opt_c1ff1.m
%   Fit the parameters of a C1-FFL to the intermediate approximation 
%
% Auxiliary file:
%   ode_ia          Use with ODE45 to solve the intermediate approximation
%   ode_C1FFL       Solves the ODE for C1-FFL
%   obj_C1FFL       Use with lsqcurvefit to fit the C1-FFL parameters 
%   plot_c1ffl_opt  Plot the graphs 
% 

% Figure 6d: response to a triangular pulse 
% Main file: main_triangular.m 
% 

