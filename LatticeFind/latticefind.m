% Function finds best lattice pars for 2Dcell with equal sides and
% LatticeAngle between them

% TODO:
% 1) extend to 3D
% 2) add varargin


function [LatPar, Angle] = latticefind(Points, LatticeAngle)

if numel(Points)>2
   error('Too many Points') 
end

PointsAngle = tan2d([Points.X], [Points.Y]);
[InitAngle,MinIndex] = min(PointsAngle);
Temp = Points(1);
Points(1) = Points(MinIndex);
Points(MinIndex) = Temp;

InitDatPar = mean(([Points.X].^2 + [Points.Y].^2).^0.5);
InitialSwarm = [InitDatPar InitAngle];

% solution finding
modelF = @(v) model(v(:,1), v(:,2), Points, LatticeAngle, false);
MaxSwarm = 20;
Options = optimoptions('particleswarm','SwarmSize',MaxSwarm, ...
    'UseParallel',true, ...
    'MaxIterations', MaxSwarm*20, ...
    'MaxStallIterations', 10, ...
    'FunctionTolerance', 1e-9, ...
    'ObjectiveLimit', -inf, ...
    'DisplayInterval', 1, ...
    'InitialSwarmMatrix', InitialSwarm,...
    'UseVectorized', false, ...
    'PlotFcn',[],'Display','none'); %PlotFcn pswplotbestf/[] Display iter/none
%         D   angle
Lower = [  0.5     0];
Upper = [ 10      360]; %FIXME create input parameter for limits
Output = particleswarm(modelF,2,Lower,Upper,Options);
clearvars Options modelF Upper Lower MaxSwarm

LatPar = Output(1);
Angle = Output(2);


    model(LatPar, Angle, Points, LatticeAngle, true);
    %disp(['LatPar = ' num2str(LatPar)])
    %disp(['Angle = ' num2str(Angle) ' deg'])


end


function out = model(LatPar, Angle, Points, LatticeAngle, DrawCmd)
px1 = Points(1).X;
py1 = Points(1).Y;
px2 = Points(2).X;
py2 = Points(2).Y;

x1 = LatPar*cosd(Angle);
y1 = LatPar*sind(Angle);
x2 = LatPar*cosd(Angle+LatticeAngle);
y2 = LatPar*sind(Angle+LatticeAngle);

Dist1 = ((x1-px1).^2 + (y1-py1).^2).^0.5;
Dist2 = ((x2-px2).^2 + (y2-py2).^2).^0.5;

out = ((Dist1).^2+(Dist2).^2);

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

function Angle = tan2d(x,y)
Angle = 90-atan2(x,y)*180/pi;
range = Angle<0;
Angle(range) = Angle(range) + 360;
end

