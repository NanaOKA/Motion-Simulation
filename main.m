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
wave_type='complex'; % Options are complex and regular
dir_deg = 80;        % wave direction [rad]

if strcmp(wave_type,'regular') == 1
    T = 9.0;             % wave period [s]
    L = 40.0;            % wavelength  [m] Must satisfy deep water condition
    A = 1.5;             % wave amplitude [m] Must satisfy steepness condition
    H = 2*A;             % wave height [m] Must satisfy steepness condition
    phi = 0;             % wave phase [rad]
    
    deep_water_condition(t_depth, L);
    
    %% SIMULATION
    dt = 0.1;                  % Time step [s], smaller = more fluid simulation
    time_sample_rate = 1/dt;   % Sampling rate [Hz]
    time_start = 0;            % Start time 
    time_end = 10;            % End time [s], length of simulation
    t = time_start:dt:time_end;

    for time_wave = time_start:dt:time_end
        %% CHECK TYPE OF SURFACE TO GENERATE
        % Create regular wave 
        H_wave = A;
        [f,w,k,c,st,dir_rad,Z,kx,ky,W,Phi] = create_wave(H_wave,L,T,dir_deg,time_wave,X,Y,phi,t_depth);

        %% CREATE FIGURE FOR WAVE AND OBJECT
        sea_surface = figure(1);
        m = mesh(X, Y, Z);

        %view([90,0])                       % default is -37.5, 30

        config_plot; 
        set(m,'FaceLighting','gouraud','FaceColor','texturemap','AmbientStrength',0.5,'EdgeColor','none');
        axis([0 t_length 0 t_width t_depth-(H+4) t_depth + (H+5)]);
        pause(0.05)
    end
    
elseif strcmp(wave_type,'complex') == 1
    spectrum='pierson'; %( Options are ....)
    Hs = 0.83;             % Significant wave height [m]
    Tp = 8.37;             % Peak/modal period [s]
    w = 0.01:0.01:5;       % Frequency [rad/s]
    
    Fp = 1/Tp;             % Peak/modal frequency [s]
    Wp = (2*pi)/Tp;        % Peak/modal frequency [rad/s]
    
    % Creat complex surface from spectrum;
    [Sw_Hs_Tp_pierson,w,Sf_Hs_Tp_pierson,f,Sw_Hs_pierson,L,T,k]=create_spectrum(w,Hs,Tp,Wp,Fp,dir_deg);

    
else
    disp('Wrong input for wave_type')
end


