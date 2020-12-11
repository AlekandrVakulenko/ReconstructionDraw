% General Test

%% General Test

% load data
load('main_3D_01.mat')

% I1D(1:2:end) = NaN;

%% function test 1
clc
% settings ------------------------
Vect = [1 1 0]; %basic vector
BasicPoint = [2 0 0]; %basic point
Trange = [-0.7 0.7]; %line parameter
Radius = 0.04;
NumberOfBins = 200;
ProjectionAxis = 'H';
%----------------------------------


[BinCoord, SumIntInBin] = CylinderFrom3D(...
    H1D, K1D, L1D, I1D, ...
    Vect, BasicPoint, Trange, ...
    Radius, NumberOfBins, ProjectionAxis);


% figure
hold on
plot(BinCoord, SumIntInBin)
set(gca,'yscale','log')
title(['Projection to ' ProjectionAxis ' axis'])

trapz(BinCoord, SumIntInBin)





%% function test 2 with debug plots
clc
% settings ------------------------
Vect = [1 0 0]; %basic vector
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

% figure
% plot(BinCoord, SumIntInBin)
% set(gca,'yscale','log')
% title(['Projection to ' ProjectionAxis ' axis'])


%%
























