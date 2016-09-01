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

[X, Y] = create_world(t_length, t_width, t_grid);

%% WAVE
wave_type='regular';

if strcmp(wave_type,'regular') == 1
    T = 9.0;             % wave period [s]
    L = 40.0;            % wavelength  [m] Must satisfy deep water condition
    A = 0.5;             % wave amplitude [m] Must satisfy steepness condition
    H = 2*A;             % wave height [m] Must satisfy steepness condition
    dir_deg = 50;        % wave direction [rad]
    phi = 0;             % wave phase [rad]
    
    deep_water_condition(t_depth, L);
elseif strcmp(wave_type,'complex') == 1
    
else
    disp('Wrong input for wave_type')
end

%% SIMULATION
dt = 0.1;                  % Time step [s], smaller = more fluid simulation
time_sample_rate = 1/dt;   % Sampling rate [Hz]
time_start = 0;            % Start time 
time_end = 140;            % End time [s], length of simulation
t = time_start:dt:time_end;

for time_wave = time_start:dt:time_end
    if strcmp(wave_type,'regular') == 1
        % Create regular wave 
        [f,w,k,c,st,dir_rad,Z,kx,ky,W,Phi] = create_wave(H_wave,L,T,dir_deg,time_wave,X,Y,phi,t_depth);
    
    elseif strcmp(wave_type,'complex') == 1
    
    else
        
    end
    
end
