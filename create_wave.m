function [f,w,k,c,st,dir_rad,Z,kx,ky,W,Phi] = create_wave(H,L,T,dir_deg,t,X,Y,phi,t_depth)
% The function takes the following inputs
% H,        wave height    [m]
% L,        wave length    [m]
% T,        wave period    [s]
% dir_deg,  wave direction [deg]
% t,        time           [s] 
% X, Y      grid 
% phi       wave phase     [rad]
% and computes and returns the following parameters
% f,        wave frequency          [Hz]
% w,        circular wave frequency [rad/s]
% k,        wave number             [rad/m]
% c,        phase velocity          [m/s] 
% st,       wave steepness          [dimensionless]  
% dir_rad,  wave direction in deg   [rad]
% Z,        surface elevation       [m]

f = 1 / T;
w = 2 * pi * f;
k = (2 * pi) / L;
c = L / T;                          % w /k
st = H/ L;
dir_rad = dir_deg * (pi / 180);
 
W = w*ones(size(X));
Phi =phi*ones(size(X));

kx = k * cosd(dir_deg);           % [rad/m]
ky = k * sind(dir_deg);           % [rad/m]

Z = H*sin((kx*X) + (ky*Y) - (W.*t) + Phi) + t_depth;
end

