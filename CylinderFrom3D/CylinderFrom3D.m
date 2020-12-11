% V1.1.0

% Release notes:
% V1.0.01:
% First release
% 
% V1.1.0:
% New faster cases for Vect == [1 0 0] or [0 1 0] or [0 0 1]
%

function [BinCoord, MeanIntInBin] = CylinderFrom3D(varargin)
narginchk(10,11);

% Matlab version control
[~, MatlabYear] = version;
MatlabYear = str2double(MatlabYear(end-3:end));


% input arguments
H1D = varargin{1}; % line array of H coord
K1D = varargin{2}; % line array of K coord
L1D = varargin{3}; % line array of L coord
I1D = varargin{4}; % line array of Intensity
Vect = varargin{5}; %(1-by-3 array) direction of cylinder
BasicPoint = varargin{6};%(1-by-3 array) basic point of cylinder T (line parameter) = 0;
Trange = varargin{7}; %(1-by-2 array) of T (line parameter) range [from to]
Radius = varargin{8}; %(scalar) radius of cylinder
NumberOfBins = varargin{9}; %(whole scalar) Number of bins in cylinder
ProjectionAxis = varargin{10}; %(single char) H, K, L, T(along cylinder)

if nargin == 11
    DrawCmd = varargin{11};%(char array) 'Draw' for debug plots (could be ignored)
else
    DrawCmd = '';
end

% ----------------------------------- arg check -----------------------------------
if find(size(Vect) ~= [1 3])
    error('Vect must be 1-by-3 vector')
end

if find(size(BasicPoint) ~= [1 3])
    error('BasicPoint must be 1-by-3 vector')
end

if find(size(Trange) ~= [1 2])
    error('Trange must be 1-by-2 vector')
end

if Radius <= 0
    error('Radius must be greater than zero')
end

if NumberOfBins <= 0
    error('NumberOfBins must be greater than zero')
end

if round(NumberOfBins) ~= NumberOfBins
    NumberOfBins = round(NumberOfBins);
    warning('NumberOfBins was rounded')
end

if ProjectionAxis~='H' && ...
        ProjectionAxis~='K' && ...
        ProjectionAxis~='L' && ...
        ProjectionAxis~='T'
    error('ProjectionAxis cloud by only ''H'', ''K'', ''L'', ''T''')
end

Nh = numel(H1D);
Nk = numel(K1D);
Nl = numel(L1D);
Ni = numel(I1D);
if Nh~=Nk || Nh~=Nl || Nh~=Ni
    error('H1D, K1D, L1D, I1D must have the same size N-by-1')
end

if size(H1D,2)~=1
    error('H1D, K1D, L1D, I1D must have the same size N-by-1')
end

if size(K1D,2)~=1
    error('H1D, K1D, L1D, I1D must have the same size N-by-1')
end

if size(L1D,2)~=1
    error('H1D, K1D, L1D, I1D must have the same size N-by-1')
end

if size(I1D,2)~=1
    error('H1D, K1D, L1D, I1D must have the same size N-by-1')
end

% ----- projection check -----
NumberOfZerosInVect = numel(find(Vect==0));
if NumberOfZerosInVect == 3
    error('Direction is zero vector')
end

Vect = Vect./norm(Vect); %input vector normalization


if ProjectionAxis == 'H'
    PhaseMargin = abs(acosd(Vect(1))-90);
    if PhaseMargin<5
        error('PhaseMargin is too bad for projection on H axis')
    end
elseif ProjectionAxis == 'K'
    PhaseMargin = abs(acosd(Vect(2))-90);
    if PhaseMargin<5
        error('PhaseMargin is too bad for projection on K axis')
    end
elseif ProjectionAxis == 'L'
    PhaseMargin = abs(acosd(Vect(3))-90);
    if PhaseMargin<5
        error('PhaseMargin is too bad for projection on L axis')
    end
end
% -----------------------------------------------------------------------



if DrawCmd == "Draw"
    Draw = 1;
    
    % Find [min max] for H, K, L
    Hmin = min(H1D);
    Hmax = max(H1D);
    Kmin = min(K1D);
    Kmax = max(K1D);
    Lmin = min(L1D);
    Lmax = max(L1D);
    
    figure
else
    Draw = 0;
end
%----------------------------------



% ------------------------------Find Cylinder------------------------------

if NumberOfZerosInVect == 2
%     'case 1'
    TransformedArray = ([H1D, K1D, L1D]' - BasicPoint');
    
    AlongIndex = find(Vect ~= 0);
    
    OrthIndex1 = mod(AlongIndex,3)+1;
    OrthIndex2 = mod(AlongIndex+1,3)+1;
    FullConversionDone = 0;
else
%     'case 2'
    % CreateOrth
    [Orth1, Orth2] = createOrth(Vect);
    Matrix = [Vect; Orth1; Orth2]; % transformation matrix / new coord = tM * vetc(col) in basic coord
    %-----------
    
    % convet coord
    TransformedArray = Matrix*([H1D, K1D, L1D]' - BasicPoint');
    % First coordinate of TransformedArray is a cylinder length (natural parameterization),
    % directed in Vect. Another two coordinates of TA are chosen arbitrarily with 90deg to main.
    %-------------
    AlongIndex = 1;
    OrthIndex1 = 2;
    OrthIndex2 = 3;
    FullConversionDone = 1;
end


% Radius condition
Condition1 = TransformedArray(OrthIndex1,:).^2 + TransformedArray(OrthIndex2,:).^2 <= Radius.^2;
% Cylinder length condition
Condition2 = TransformedArray(AlongIndex,:)>Trange(1) & TransformedArray(AlongIndex,:)<Trange(2);
% full condition
Condition = Condition1 & Condition2;
%-----------------

% TransformedArray rebuild by condition
TransformedArray = TransformedArray(:,Condition);
% Intensity in cylinder
PartIntensity = I1D(Condition);

% Draw cylinder
if Draw == 1
    subplot(1,2,1)
    Partxg = H1D(Condition);
    Partyg = K1D(Condition);
    Partzg = L1D(Condition);
    plot3(Partxg, Partyg, Partzg,'.')
    xlabel('H')
    ylabel('K')
    zlabel('L')
    xlim([Hmin Hmax])
    ylim([Kmin Kmax])
    zlim([Lmin Lmax])
    title('Whole cylinder')
    drawnow
    % axis equal
end

% -------------------------place Bins in Cylinder--------------------------
BinEdges = linspace(Trange(1),Trange(2),NumberOfBins+1);
BinCenters = (BinEdges(1:end-1)+BinEdges(2:end))/2;

if Draw == 1
    subplot(1,2,2)
    hold on
end

MeanIntInBin = zeros(1,NumberOfBins);
for index = 1:NumberOfBins
    
    %current BinEdges
    BinStart = BinEdges(index);
    BinEnd = BinEdges(index+1);
    
    % Bin condition (Сylinder length from BinStart to BinEnd)
    Condition = TransformedArray(AlongIndex,:)>=BinStart & TransformedArray(AlongIndex,:)<BinEnd;
    IntensityInBin = PartIntensity(Condition);
    
    NaNcondition = isnan(IntensityInBin);
    if MatlabYear < 2019
        MeanIntInBin(index) = mean(IntensityInBin(~NaNcondition));
    elseif MatlabYear >= 2019
        MeanIntInBin(index) = mean(IntensityInBin(~NaNcondition),'all');
    end
    
    % Draw Bins
    if Draw == 1
        HinBin = Partxg(Condition);
        KinBin = Partyg(Condition);
        LinBin = Partzg(Condition);
        plot3(HinBin, KinBin, LinBin,'.')
        xlabel('H')
        ylabel('K')
        zlabel('L')
        xlim([Hmin Hmax])
        ylim([Kmin Kmax])
        zlim([Lmin Lmax])
        title('Bins of cylinder')
        drawnow
    end
end
% -------------------------------------------------------------------------

% zeroing negative
MeanIntInBin(MeanIntInBin<=0) = 0;

% find NaN in output
NaNcondition = isnan(MeanIntInBin);
BinCenters(NaNcondition) = [];
MeanIntInBin(NaNcondition) = [];


% --------------- find projections of BinCenters to H, K, L ---------------
if FullConversionDone == 1
    MidValue = ([BinCenters; zeros(size(BinCenters)); zeros(size(BinCenters))])'/Matrix' + BasicPoint;
    OutHgrid = MidValue(:,1)';
    OutKgrid = MidValue(:,2)';
    OutLgrid = MidValue(:,3)';
else
    MidMatrix = zeros(3,numel(BinCenters));
    MidMatrix(AlongIndex,:) = BinCenters;
    
    MidValue = (MidMatrix)' + BasicPoint;
    OutHgrid = MidValue(:,1)';
    OutKgrid = MidValue(:,2)';
    OutLgrid = MidValue(:,3)';
end
% -------------------------------------------------------------------------


% Draw BinCenters points
if Draw == 1
    plot3(OutHgrid,OutKgrid,OutLgrid,'o','markersize',8,'markerfacecolor','b')
end


% Tstep,Hstep,Kstep,Lstep are steps between points along T,H,K,L
Tstep = mean(abs(diff(BinCenters)));
Hstep = mean(abs(diff(OutHgrid)));
Kstep = mean(abs(diff(OutKgrid)));
Lstep = mean(abs(diff(OutLgrid)));

if ProjectionAxis == 'H'
    BinCoord = OutHgrid; %output grid
    MeanIntInBin = MeanIntInBin*(pi*Radius^2)*Tstep/Hstep; %integration
    
elseif ProjectionAxis == 'K'
    BinCoord = OutKgrid; %output grid
    MeanIntInBin = MeanIntInBin*(pi*Radius^2)*Tstep/Kstep; %integration
    
elseif ProjectionAxis == 'L'
    BinCoord = OutLgrid; %output grid
    MeanIntInBin = MeanIntInBin*(pi*Radius^2)*Tstep/Lstep; %integration
    
elseif ProjectionAxis == 'T'
    BinCoord = BinCenters; %output grid
    MeanIntInBin = MeanIntInBin*(pi*Radius^2); %integration
end


end








