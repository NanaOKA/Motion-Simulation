function [] = config_plot()
% Additional configurations to make plot look realistic
% view([0,0])                         % default is -37.5, 30
set(gca,'DataAspectRatio',[1 1 1])  % default [3 1 1]
xlabel('X, Length axis ');
ylabel('Y, Width axis');
zlabel('Z, Height axis');

%set(gcf, 'Position',[0 0 1000 1000]);

light('Position',[0,0,20],'Style','infinite')
colormap winter
%shading interp (CANNOT BE USED IN COMBINATION WITH PATCH)
zlim('manual')
%pause(0.05)
end

