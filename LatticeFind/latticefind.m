% Function finds best lattice pars for cubic cell on 2D

% TODO:
% 1) extend to 3D
% 2) add varargin


function [LatPar, Angle] = latticefind(Points, DrawCmd)

% solution finding
modelF = @(v) model(v(:,1), v(:,2), Points, false);
MaxSwarm = 10;
Options = optimoptions('particleswarm','SwarmSize',MaxSwarm, ...
    'UseParallel',true, ...
    'MaxIterations', MaxSwarm*20, ...
    'MaxStallIterations', 10, ...
    'FunctionTolerance', 1e-6, ...
    'ObjectiveLimit', -inf, ...
    'DisplayInterval', 1, ...
    'UseVectorized', true, ...
    'PlotFcn',[],'Display','none'); %PlotFcn pswplotbestf/[] Display iter/none
%         D   angle
Lower = [  0.01     0];
Upper = [ 10      360]; %FIXME create input parameter for limits
Output = particleswarm(modelF,2,Lower,Upper,Options);
clearvars Options modelF Upper Lower MaxSwarm

LatPar = Output(1);
Angle = Output(2);

if DrawCmd
    model(LatPar, Angle, Points, true);
    disp(['LatPar = ' num2str(LatPar)])
    disp(['Angle = ' num2str(Angle) ' deg'])
end

end


function out = model(LatPar, Angle, Points, DrawCmd)
px1 = Points(1).X;
py1 = Points(1).Y;
px2 = Points(2).X;
py2 = Points(2).Y;

x1 = LatPar.*cosd(Angle);
y1 = LatPar.*sind(Angle);
x2 = LatPar.*cosd(Angle+90);
y2 = LatPar.*sind(Angle+90);

Dist1 = ((x1-px1).^2 + (y1-py1).^2).^0.5;
Dist2 = ((x2-px2).^2 + (y2-py2).^2).^0.5;

W1 = 10;
W2 = 1;

out = ( W1*(Dist1 - Dist2).^2 + W2*((Dist1).^2+(Dist1).^2) );

if DrawCmd
    figure
    hold on
    
    plot([x1 0], [y1 0], 'b' ,'linewidth', 1)
    plot([x2 0], [y2 0], 'b' ,'linewidth', 1)
    plot(px1, py1, 'r.', 'markersize', 12)
    plot(px2, py2, 'r.', 'markersize', 12)
    xlim([-3 3])
    ylim([-3 3])
    xline(0);
    yline(0);
    axis equal
    
end
end
