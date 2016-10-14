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
f_actual = [ 0.        ,  0.00853333,  0.01706667,  0.0256    ,  0.03413333, ...
        0.04266667,  0.0512    ,  0.05973333,  0.06826667,  0.0768    , ...
        0.08533333,  0.09386667,  0.1024    ,  0.11093333,  0.11946667, ...
        0.128     ,  0.13653333,  0.14506667,  0.1536    ,  0.16213333, ...
        0.17066667,  0.1792    ,  0.18773333,  0.19626667,  0.2048    , ...
        0.21333333,  0.22186667,  0.2304    ,  0.23893333,  0.24746667, ...
        0.256     ,  0.26453333,  0.27306667,  0.2816    ,  0.29013333, ...
        0.29866667,  0.3072    ,  0.31573333,  0.32426667,  0.3328    , ...
        0.34133333,  0.34986667,  0.3584    ,  0.36693333,  0.37546667, ...
        0.384     ,  0.39253333,  0.40106667,  0.4096    ,  0.41813333, ...
        0.42666667,  0.4352    ,  0.44373333,  0.45226667,  0.4608    , ...
        0.46933333,  0.47786667,  0.4864    ,  0.49493333,  0.50346667, ...
        0.512     ,  0.52053333,  0.52906667,  0.5376    ,  0.54613333, ...
        0.55466667,  0.5632    ,  0.57173333,  0.58026667,  0.5888    , ...
        0.59733333,  0.60586667,  0.6144    ,  0.62293333,  0.63146667, ...
        0.64      ];
    
Sf_actual = [  4.80933967e-04,   2.59968959e-04,   1.80607171e-05, ...
         8.23803017e-05,   9.27381702e-04,   1.97444832e-03, ...
         4.98651156e-03,   1.94989991e-02,   2.93637458e-02, ...
         3.00313192e-02,   4.93691502e-02,   1.35301360e-01, ...
         2.72210594e-01,   3.99030786e-01,   5.87927513e-01, ...
         4.87578858e-01,   3.36700935e-01,   3.23199934e-01, ...
         2.63145331e-01,   2.36153224e-01,   2.39314388e-01, ...
         2.62295005e-01,   2.07925857e-01,   2.05850112e-01, ...
         1.62826700e-01,   1.24494574e-01,   8.37922010e-02, ...
         9.62941979e-02,   5.86984694e-02,   5.73607467e-02, ...
         6.22845540e-02,   5.02356024e-02,   3.07947839e-02, ...
         2.26403507e-02,   2.49004958e-02,   2.16043504e-02, ...
         1.85200587e-02,   1.45731238e-02,   1.66720601e-02, ... 
         1.36901130e-02,   1.17402402e-02,   5.85627617e-03, ...
         7.65468065e-03,   8.96966160e-03,   6.81397361e-03, ...
         5.76646795e-03,   6.94438103e-03,   6.52383065e-03, ...
         5.13610471e-03,   4.56267389e-03,   3.78358954e-03, ...
         3.16105134e-03,   2.60600111e-03,   2.75178829e-03, ...
         2.36757082e-03,   2.17832810e-03,   2.24448600e-03, ...
         2.50232400e-03,   2.35799631e-03,   2.01414438e-03, ...
         1.72372318e-03,   1.65979449e-03,   2.15241928e-03, ...
         2.48461202e-03,   2.09848616e-03,   2.30355194e-03, ...
         2.44817751e-03,   2.68787299e-03,   2.60007992e-03, ...
         2.08554053e-03,   1.65139806e-03,   1.55937598e-03, ...
         4.80346577e-04,   2.05761750e-05,   9.69270053e-06, ... 
         5.61453674e-06];

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

Sf_Hs_Tp_pierson = 2*pi*af.*bf.*cf.*exp(df.*ef);
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

Sf_Hs_Tp_bret = 2*pi*B1f.*B2f.*(Hs).*exp(B3f.*B4f);


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

m0_actual = trapz((f_actual.^0).*Sf_actual);
m1_actual = trapz((f_actual.^1).*Sf_actual);
m2_actual = trapz((f_actual.^2).*Sf_actual);
m4_actual = trapz((f_actual.^4).*Sf_actual);

[M_pier,I_pier]=max(Sf_Hs_Tp_pierson);
[M_bret,I_bret]=max(Sf_Hs_Tp_bret);
[M_actual,I_actual]=max(Sf_actual);

Calc_Tp_pier = 1/f(I_pier)     % Peak period
Calc_Tp_bret = 1/f(I_bret)
Calc_Tp_actual = 1/f_actual(I_actual)
 
Calc_Hs_pier = 0.4*sqrt(m0_pier/(2*pi))   % Significant wave height [m]
Calc_Hs_bret = 0.4*sqrt(m0_bret/(2*pi))   % Division by 2 is to do with conversion from w to f
Calc_Hs_actual = 0.4*sqrt(m0_actual)

Calc_Ta_pier = m0_pier/m1_pier          % Average wave period [s] = Mean centroid wave
Calc_Ta_bret = m0_bret/m1_bret 
Calc_Ta_actual = m0_actual/m1_actual

Calc_Tz_pier = sqrt(m0_pier/m2_pier)    % Mean zero-crossing wave period = Significant wave period, also known as Ts
Calc_Tz_bret = sqrt(m0_bret/m2_bret)
Calc_Tz_actual = sqrt(m0_actual/m2_actual)

Calc_T_p_m_pier = sqrt(m2_pier/m4_pier) % Mean period between maxima 
Calc_T_p_m_bret = sqrt(m2_bret/m4_bret)
Calc_T_p_m_actual = sqrt(m2_actual/m4_actual)

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

% % See http://www.ultramarine.com/hdesk/document/papers/sea_spectra_simplified.pdf
% % INPUTS FOR WAVE HEIGHTS IN FEET
% figure;
% plot(w,Sw_Hs_Tz_bret);
% xlabel('Frequency, \omega [rad/s]');
% ylabel('Height Double Spectrum [ft^2s]');
% title('Bretschneider Spectrum(wrt rad/s) taking Hs and Tz as inputs')
% legend('BRET (Hs, Tz)')
% %hold on;

figure;
plot(w,Sw_Hs_Tp_pierson,'r',w,Sw_Hs_Tp_bret,'b',w,Sw_Hs_pierson,'r--');
xlabel('Frequency, \omega [rad/s]');
ylabel('Power Spectral Density [m^2s/rad]');
legend('Pierson-Moskowitz','Bretschneider','PM (only Hs)');


figure;
plot(f,(Sf_Hs_Tp_pierson),'r',f,(Sf_Hs_Tp_bret),'b');
xlabel('Frequency, f [Hz]');
ylabel('Power Spectral Density [m^2/Hz]');
legend('Pierson-Moskowitz','Bretschneider');
hold on

plot(f_actual,Sf_actual,'--');

end


