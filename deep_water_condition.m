function [] = deep_water_condition(t_depth, L)
% Checks whether the deep water condition is satisfied
if t_depth < (L/2)
    disp(fprintf(['Deep water condition not satisfied.\n', ...
        'Depth of water is less than half the wavelength.\n' ...
        'Either increase water depth or decrease wavelength.']))
    return
end

end

