function [Sw,w,Sf,f] = create_spectrum(Hs,Tp,dir_deg)
% This function takes the following inputs
% Hs,       significant wave height[m]
% Tp,       peak/modal wave period[s]
% and returns
% Sw,       Power spectral density
% w,        Frequency [rad/s]
% f,        Frequency [Hz]

w = 0.01:0.01:5;  % Frequency [rad/s]
wx = w .* cosd(dir_deg);
wy = w .* sind(dir_deg);

aw = 5*(pi^4);
bw = (Hs^2)/(Tp^4);
cw = (sqrt((wx.^2)+(wy.^2)).^-5);
dw = (-20*(pi^4))/(Tp^4);
ew = (sqrt((wx.^2)+(wy.^2)).^-4);

Sw = aw.*bw.*cw.*exp(dw.*ew);

f = w./(2*pi); % Frequency [Hz]
fx = f .* cosd(dir_deg);
fy = f .* sind(dir_deg);

af = 5/(32*pi);
bf = (Hs^2)/(Tp^4);
cf = (sqrt((fx.^2)+(fy.^2)).^-5);
df = (-20/(16*(Tp^4)));
ef = (sqrt((fx.^2)+(fy.^2)).^-4);

Sf = af.*bf.*cf.*exp(df.*ef);

m0 = trapz((f.^0).*Sf);
m1 = trapz((f.^1).*Sw);
m2 = trapz((f.^2).*Sw);
m4 = trapz((f.^4).*Sw);

Calc_Hs = 0.4*sqrt(m0);   % Significant wave height [m]
Calc_Ta = m0/m1;          % Average wave period [s] - Mean centroid wave
Calc_Tz = sqrt(m0/m2);    % Mean zero-crossing wave period
Calc_T_p_m = sqrt(m2/m4) % Mean period between maxima 

figure;
plot(w,Sw);
xlabel('Frequency, \omega [rad/s]');
ylabel('Power Spectral Density [m^2s/2 \pi rad]');

figure;
plot(f,Sf);
xlabel('Frequency, f [Hz]');
ylabel('Power Spectral Density [m^2/Hz]');
end


