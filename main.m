%%%%%%%%%%%%%%%%%%%%%
% Numerical Time Domain Simulation
% to investigate motions of a fishing vessel
%
% Nana Okra Abankwa
% n.abankwa@soton.ac.uk
% September 2016
%%%%%%%%%%%%%%%%%%%%%

clc       % clear all input and output from command window display
clear     % clear all variables from workspace
close all % close all open figure windows

%% MISCELLANEOUS USER INPUT
w_density = 1029;   % [kg/m^3] Density of water
g = 9.81;           % [Newtons/kg] Gravitational force

%% WORLD
t_length = 30.0;    % Tank length [m]  default = 100 X axis
t_width  = 30.0;    % Tank width [m]  default = 100 Y axis
t_depth  = 30.0;    % Tank length [m] default = 100 Z axis 
t_grid   = 0.05;    % Size of grid tank will be divided into [m]