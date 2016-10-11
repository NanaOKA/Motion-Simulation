function [Sw_Hs_Tp_pierson,w,Sf_Hs_Tp_pierson,f,Sw_Hs_pierson,L,T,k] = create_spectrum(w,Hs,Tp,Wp,dir_deg)
% This function takes the following inputs
% Hs,       significant wave height[m]
% Tp,       peak/modal wave period[s]
% and returns
% Sw,       Power spectral density
% w,        Frequency [rad/s]
% f,        Frequency [Hz]

g = 9.80665;         % Acceleration due to gravity [m/s^2]
L = 2*pi*g*(w.^-2);  % Wavelength [L]
T = 2*pi./(w);       % Wave period[s]
k = (2*pi)./L;       % Wave number [rad/m]


%% PIERSON-MOSKOWITZ SPECTRUM
% See http://www.codecogs.com/library/engineering/fluid_mechanics/waves/spectra/pierson_moskowitz.php
wx = w .* cosd(dir_deg);
wy = w .* sind(dir_deg);

aw = 5*(pi^4);
bw = (Hs^2)/(Tp^4);
cw = (w.^-5);
dw = (-20*(pi^4))/(Tp^4);
ew = (w.^-4);

Sw_Hs_Tp_pierson = aw.*bw.*cw.*exp(dw.*ew);

f = w./(2*pi); % Frequency [Hz]
fx = f .* cosd(dir_deg);
fy = f .* sind(dir_deg);

af = 5/(32*pi);
bf = (Hs^2)/(Tp^4);
cf = (f.^-5);
df = (-20/(16*(Tp^4)));
ef = (f.^-4);

Sf_Hs_Tp_pierson = af.*bf.*cf.*exp(df.*ef);

%See https://ocw.mit.edu/courses/mechanical-engineering/2-22-design-principles-for-ocean-vehicles-13-42-spring-2005/readings/r8_wavespectra.pdf
a1 = 0.0081;
b1 = (g^2)*(w.^-5);
c1 = -0.032;
d1 = g./((w.^2)*Hs);

Sw_Hs_pierson = a1.*b1.*exp(c1*d1.^2);
%S1_pierson same as Sw_pierson with Tp~4.2, Hs~0.83

%% BRETSCHNEIDER SPECTRUM (ITTC Standard)
B1 = 1.25/4;
B2 = (Wp^4)./(w.^5);
B3 = -1.25;
B4 = (Wp./w).^4;

Sw_Hs_Tp_bret = B1.*B2.*Hs.*exp(B3.*B4);
%% SPECTRAL MOMENTS
m0 = trapz((f.^0).*Sf_Hs_Tp_pierson);
m1 = trapz((f.^1).*Sf_Hs_Tp_pierson);
m2 = trapz((f.^2).*Sf_Hs_Tp_pierson);
m4 = trapz((f.^4).*Sf_Hs_Tp_pierson);

Calc_Hs = 0.4*sqrt(m0);   % Significant wave height [m]
Calc_Ta = m0/m1;          % Average wave period [s] - Mean centroid wave
Calc_Tz = sqrt(m0/m2);    % Mean zero-crossing wave period
Calc_T_p_m = sqrt(m2/m4); % Mean period between maxima 

figure;
plot(f,Sf_Hs_Tp_pierson);
xlabel('Frequency, f [Hz]');
ylabel('Power Spectral Density [m^2/Hz]');
title('Pierson-Moskowitz Spectrum(wrt Hz) taking Hs and Tp as inputs')
%hold on;

figure;
plot(w,Sw_Hs_Tp_pierson);
xlabel('Frequency, \omega [rad/s]');
ylabel('Power Spectral Density [m^2s/2 \pi rad]');
title('Pierson-Moskowitz Spectrum(wrt rad/s) taking Hs and Tp as inputs')
legend('PM (Hs,Tp)');
hold on;

figure;
plot(w,Sw_Hs_pierson);
xlabel('Frequency, \omega [rad/s]');
ylabel('Power Spectral Density [m^2s/2 \pi rad]');
title('Pierson-Moskowitz Spectrum(wrt rad/s) taking taking only Hs as an inputs. Tp incorprated in equation')
%hold on;

figure;
plot(w,Sw_Hs_Tp_bret);
xlabel('Frequency, \omega [rad/s]');
ylabel('Power Spectral Density [m^2s/2 \pi rad]');
title('Bretschneider Spectrum(wrt rad/s) taking Hs and Tp as inputs')
legend('BRET (Hs, Tp)')
%hold on;
end


