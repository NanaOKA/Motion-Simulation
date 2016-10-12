function [Sw_Hs_Tp_pierson,w,Sf_Hs_Tp_pierson,f,Sw_Hs_pierson,L,T,k] = create_spectrum(w,Hs,Tp,Wp,Fp,dir_deg)
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
% 
% %See https://ocw.mit.edu/courses/mechanical-engineering/2-22-design-principles-for-ocean-vehicles-13-42-spring-2005/readings/r8_wavespectra.pdf
a1 = 0.0081;
b1 = (g^2)*(w.^-5);
c1 = -0.032;
d1 = g./((w.^2)*Hs);

Sw_Hs_pierson = a1.*b1.*exp(c1*d1.^2);
% %S1_pierson same as Sw_pierson with Tp~4.2, Hs~0.83

%% BRETSCHNEIDER SPECTRUM (ITTC Standard)
% See https://ocw.mit.edu/courses/mechanical-engineering/2-22-design-principles-for-ocean-vehicles-13-42-spring-2005/readings/r8_wavespectra.pdf
B1 = 1.25/4;
B2 = (Wp^4)./(w.^5);
B3 = -1.25;
B4 = (Wp./w).^4;

Sw_Hs_Tp_bret = B1.*B2.*Hs.*exp(B3.*B4);

B1f = 1.25/(4*2*pi);
B2f = (Fp^4)./(f.^5);
B3f = -1.25;
B4f = (Fp./f).^4;

Sf_Hs_Tp_bret = B1f.*B2f.*Hs.*exp(B3f.*B4f);


% % See http://www.ultramarine.com/hdesk/document/papers/sea_spectra_simplified.pdf
% % INPUTS FOR WAVE HEIGHTS IN FEET
% Tz = 1;
% B1_ = 4200*(Hs^2);
% B2_ = (Tz^4)*(w.^5);
% B3_ = -1050;
% B4_ = (Tz^4).*(w.^4);
% Sw_Hs_Tz_bret = B1_./(B2_.*exp(B3_./B4_));

%% SPECTRAL MOMENTS
m0_pier = trapz((f.^0).*Sf_Hs_Tp_pierson);
m1_pier = trapz((f.^1).*Sf_Hs_Tp_pierson);
m2_pier = trapz((f.^2).*Sf_Hs_Tp_pierson);
m4_pier = trapz((f.^4).*Sf_Hs_Tp_pierson);

m0_bret = trapz((f.^0).*Sf_Hs_Tp_bret);
m1_bret = trapz((f.^1).*Sf_Hs_Tp_bret);
m2_bret = trapz((f.^2).*Sf_Hs_Tp_bret);
m4_bret = trapz((f.^4).*Sf_Hs_Tp_bret);

Calc_Hs_pier = 0.4*sqrt(m0_pier)   % Significant wave height [m]
Calc_Hs_bret = 0.4*sqrt(m0_bret)

Calc_Ta_pier = m0_pier/m1_pier          % Average wave period [s] = Mean centroid wave
Calc_Ta_bret = m0_bret/m1_bret 

Calc_Tz_pier = sqrt(m0_pier/m2_pier)    % Mean zero-crossing wave period = Significant wave period, also known as Ts
Calc_Tz_bret = sqrt(m0_bret/m2_bret)

Calc_T_p_m_pier = sqrt(m2_pier/m4_pier) % Mean period between maxima 
Calc_T_p_m_bret = sqrt(m2_bret/m4_bret)


figure;
plot(f,Sf_Hs_Tp_pierson,'r:');
xlabel('Frequency, f [Hz]');
ylabel('Power Spectral Density [m^2/Hz]');
title('Pierson-Moskowitz Spectrum(wrt Hz) taking Hs and Tp as inputs')
%hold on;

figure;
plot(w,Sw_Hs_Tp_pierson,'r');
xlabel('Frequency, \omega [rad/s]');
ylabel('Power Spectral Density [m^2s/2 \pi rad]');
title('Pierson-Moskowitz Spectrum(wrt rad/s) taking Hs and Tp as inputs')
legend('PM (Hs,Tp)');
%hold on;

figure;
plot(w,Sw_Hs_pierson,'r--');
xlabel('Frequency, \omega [rad/s]');
ylabel('Power Spectral Density [m^2s/2 \pi rad]');
title('Pierson-Moskowitz Spectrum(wrt rad/s) taking taking only Hs as an inputs. Tp incorprated in equation')
%hold on;

figure;
plot(w,Sw_Hs_Tp_bret,'b');
xlabel('Frequency, \omega [rad/s]');
ylabel('Power Spectral Density [m^2s/2 \pi rad]');
title('Bretschneider Spectrum(wrt rad/s) taking Hs and Tp as inputs')
legend('BRET (Hs, Tp)')

figure;
plot(f,Sf_Hs_Tp_bret,'b:');
xlabel('Frequency, f [Hz]');
ylabel('Power Spectral Density [m^2/Hz]');
title('Bretschneider Spectrum(wrt Hz) taking Hs and Tp as inputs')
%hold on;

% figure;
% plot(w,Sw_Hs_Tz_bret);
% xlabel('Frequency, \omega [rad/s]');
% ylabel('Height Double Spectrum [ft^2s]');
% title('Bretschneider Spectrum(wrt rad/s) taking Hs and Tz as inputs')
% legend('BRET (Hs, Tz)')
% %hold on;
end


