function [Sw,w] = create_spectrum(Hs,Tp,dir_deg)
% This function takes the following inputs
% Hs,       significant wave height[m]
% Tp,       peak/modal wave period[s]
% and returns
% Sw,       Power spectral density
% w,        Frequency [rad/s]

w = 0:0.01:5;
f = w./(2*pi);

wx = w .* cosd(dir_deg);
wy = w .* sind(dir_deg);

a = 5*(pi^4);
b = (Hs^2)/(Tp^4);
c = (w.^-5);
d = (-20*(pi^4))/(Tp^4);
e = (sqrt((wx.^2)+(wy.^2)).^-5);

Sw = a.*b.*c.*exp(d.*e);
figure;
plot(w,Sw);
xlabel('Frequency, \omega [rad/s]');
ylabel('Power Spectral Density [m^2/Hz]');

figure;
plot(f,Sw);
xlabel('Frequency, f [Hz]');
ylabel('Power Spectral Density [m^2s/2 \pi rad]');
end


