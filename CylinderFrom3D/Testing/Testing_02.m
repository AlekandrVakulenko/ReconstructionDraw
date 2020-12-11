% Test for integration

%%


N = 200;

Hrange = linspace(-0.1, 1.1, N);
Krange = linspace(-0.1, 1.1, N);
Lrange = linspace(-0.6, 0.6, N);


[Hgrid, Kgrid, Lgrid] = meshgrid(Hrange, Krange, Lrange);



Peak = @(x,y,z,x0,y0,z0,Amp,Size) Amp/(pi*Size)^(3*0.5)*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)/Size);


Peak_1_H = 0.5;
Peak_1_K = 0.5;
Peak_1_L = 0.0;
Peak_1_Amp = 1;
Peak_1_Size = 0.0001;

Peak_2_H = 0.25;
Peak_2_K = 0.25;
Peak_2_L = 0.0;
Peak_2_Amp = 2;
Peak_2_Size = 0.0001;


Iall = Peak(Hgrid, Kgrid, Lgrid, Peak_1_H, Peak_1_K, Peak_1_L, Peak_1_Amp, Peak_1_Size) + ...
       Peak(Hgrid, Kgrid, Lgrid, Peak_2_H, Peak_2_K, Peak_2_L, Peak_2_Amp, Peak_2_Size);


isosurface(Hgrid, Kgrid, Lgrid, Iall, 0.05)

axis equal
xlim([-0.1, 1.1])
ylim([-0.1, 1.1])
zlim([-0.6, 0.6])
xlabel('H')
ylabel('K')
zlabel('L')



I1 = trapz(Hrange, Iall);
I2 = trapz(Krange, I1);
I3 = trapz(Lrange, I2);

disp(['full integral = ' num2str(I3) ' (2+1)'])





H1D = reshape(Hgrid, numel(Hgrid), 1);
K1D = reshape(Kgrid, numel(Kgrid), 1);
L1D = reshape(Lgrid, numel(Lgrid), 1);
I1D = reshape(Iall, numel(Iall), 1);



%% Integral of (0.5 0.5 0)
clc
% settings ------------------------
Vect = [1 0.001 0]; %basic vector
BasicPoint = [0.5 0.5 0]; %basic point
Trange = [-1.7 1.7]; %line parameter
Radius = 0.08;
NumberOfBins = 500;
ProjectionAxis = 'H';
%----------------------------------


[BinCoord, SumIntInBin] = CylinderFrom3D(...
    H1D, K1D, L1D, I1D, ...
    Vect, BasicPoint, Trange, ...
    Radius, NumberOfBins, ProjectionAxis);


% figure
hold on
plot(BinCoord, SumIntInBin)
% set(gca,'yscale','log')
title(['Projection to ' ProjectionAxis ' axis'])

trapz(BinCoord, SumIntInBin)


%% Integral of (0.5 0.5 0) + (0.25 0.25 0)
clc
% settings ------------------------
Vect = [1 1 0]; %basic vector
BasicPoint = [0.5 0.5 0]; %basic point
Trange = [-1.7 1.7]; %line parameter
Radius = 0.08;
NumberOfBins = 500;
ProjectionAxis = 'H';
%----------------------------------


[BinCoord, SumIntInBin] = CylinderFrom3D(...
    H1D, K1D, L1D, I1D, ...
    Vect, BasicPoint, Trange, ...
    Radius, NumberOfBins, ProjectionAxis);


% figure
hold on
plot(BinCoord, SumIntInBin)
% set(gca,'yscale','log')
title(['Projection to ' ProjectionAxis ' axis'])

trapz(BinCoord, SumIntInBin)










































