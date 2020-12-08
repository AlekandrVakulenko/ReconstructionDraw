

%%

% load data
load('main_3D_01.mat')


%% function test 1
clc
% settings ------------------------
Vect = [-2 2 0]; %basic vector
BasicPoint = [2 0 0]; %basic point
Trange = [-0.3 0.3]; %line parameter
Radius = 0.04;
NumberOfBins = 20;
ProjectionAxis = 'H';
%----------------------------------


[BinCoord, SumIntInBin] = CylinderFrom3D(...
    H1D, K1D, L1D, I1D, ...
    Vect, BasicPoint, Trange, ...
    Radius, NumberOfBins, ProjectionAxis);


figure
plot(BinCoord, SumIntInBin)
set(gca,'yscale','log')
title(['Projection to ' ProjectionAxis ' axis'])




%% function test 2 with debug plots
clc
% settings ------------------------
Vect = [-2 2 0]; %basic vector
BasicPoint = [2 0 0]; %basic point
Trange = [-0.3 0.3]; %line parameter
Radius = 0.04;
NumberOfBins = 20;
ProjectionAxis = 'H';
DrawCmd = 'Draw';
%----------------------------------

% now DrawCmd = 'Draw' passes as 11th function argument

[BinCoord, SumIntInBin] = CylinderFrom3D(...
    H1D, K1D, L1D, I1D, ...
    Vect, BasicPoint, Trange, ...
    Radius, NumberOfBins, ProjectionAxis, ...
    DrawCmd);


figure
plot(BinCoord, SumIntInBin)
set(gca,'yscale','log')
title(['Projection to ' ProjectionAxis ' axis'])


%%

























