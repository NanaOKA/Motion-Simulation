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
wave_type='complex';

if strcmp(wave_type,'regular') == 1
    T = 9.0;             % wave period [s]
    L = 40.0;            % wavelength  [m] Must satisfy deep water condition
    A = 1.5;             % wave amplitude [m] Must satisfy steepness condition
    H = 2*A;             % wave height [m] Must satisfy steepness condition
    dir_deg = 40;        % wave direction [rad]
    phi = 0;             % wave phase [rad]
    
    deep_water_condition(t_depth, L);
    
elseif strcmp(wave_type,'complex') == 1
    spectrum='pierson'; %( Options are ....)
    Hs = 4;             % Significant wave height
    Tp = 8;             % Peak period
    dir_deg = 10;        % wave direction [rad]
    
else
    disp('Wrong input for wave_type')
end

%% SIMULATION
dt = 0.1;                  % Time step [s], smaller = more fluid simulation
time_sample_rate = 1/dt;   % Sampling rate [Hz]
time_start = 0;            % Start time 
time_end = 10;            % End time [s], length of simulation
t = time_start:dt:time_end;

for time_wave = time_start:dt:time_end
    %% CHECK TYPE OF SURFACE TO GENERATE
    if strcmp(wave_type,'regular') == 1
        % Create regular wave 
        H_wave = A;
        [f,w,k,c,st,dir_rad,Z,kx,ky,W,Phi] = create_wave(H_wave,L,T,dir_deg,time_wave,X,Y,phi,t_depth);
    
    elseif strcmp(wave_type,'complex') == 1
        % Creat complex surface from spectrum;
        % create_spectrum()
        w = 0:0.01:5;
        
        wx = w .* cosd(dir_deg);
        wy = w .* sind(dir_deg);
        
        a = 5*(pi^4);
        b = (Hs^2)/(Tp^4);
        c = (w.^-5);
        d = (-20*(pi^4))/(Tp^4);
        e = (w.^-5);
        e = (sqrt((wx.^2)+(wy.^2)).^-5);

        Sw = a.*b.*c.*exp(d.*e);

        plot(w,Sw);
    
    else
        
    end
    
    %% CREATE FIGURE FOR WAVE AND OBJECT
    sea_surface = figure(1);
    m = mesh(X, Y, Z);
  
    %view([90,0])                       % default is -37.5, 30
    
    config_plot; 
    set(m,'FaceLighting','gouraud','FaceColor','texturemap','AmbientStrength',0.5,'EdgeColor','none');
    axis([0 t_length 0 t_width t_depth-(H+4) t_depth + (H+5)]);
    pause(0.05)
end
