
function [BinCoord, SumIntInBin] = CylinderFrom3D(varargin)
narginchk(10,11);

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

if ProjectionAxis~='H' & ...
        ProjectionAxis~='K' & ...
        ProjectionAxis~='L' & ...
        ProjectionAxis~='T'
    error('ProjectionAxis cloud by only ''H'', ''K'', ''L'', ''T''')
end

Nh = numel(H1D);
Nk = numel(K1D);
Nl = numel(L1D);
Ni = numel(I1D);
if Nh~=Nk | Nh~=Nl | Nh~=Ni
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
% -----------------------------------------------------------------------



if DrawCmd == "Draw"
    Draw = 1;
    figure
else
    Draw = 0;
end


% Find [min max] for H, K, L
Hmin = min(H1D);
Hmax = max(H1D);
Kmin = min(K1D);
Kmax = max(K1D);
Lmin = min(L1D);
Lmax = max(L1D);
% DEBUG 3 lines
[Hmin Hmax];
[Kmin Kmax];
[Lmin Lmax];
%----------------------------------


% ------------------------------Find Cylinder------------------------------
% CreateOrth
Vect = Vect./norm(Vect);
[orth1 orth2] = createOrth(Vect);
Matrix = [Vect; orth1; orth2]; % transformation matrix / new coord = tM * vetc(col) in basic coord
%-----------

% convet coord
TransformedArray = Matrix*([H1D, K1D, L1D]' - BasicPoint');
% First coordinate of TransformedArray is a cylinder length (natural parameterization),
% directed in Vect. Another two coordinates of TA are chosen arbitrarily with 90deg to main.
%-------------

% Radius condition
Condition1 = TransformedArray(2,:).^2 + TransformedArray(3,:).^2 <= Radius.^2;
% Сylinder length condition
Condition2 = TransformedArray(1,:)>Trange(1) & TransformedArray(1,:)<Trange(2);
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

SumIntInBin = [];
for index = 1:NumberOfBins
    
    %current BinEdges
    BinStart = BinEdges(index);
    BinEnd = BinEdges(index+1);
    
    % Bin condition (Сylinder length from BinStart to BinEnd)
    Condition = TransformedArray(1,:)>BinStart & TransformedArray(1,:)<BinEnd;
    
    IntensityInBin = PartIntensity(Condition);
    SumIntInBin(index) = sum(IntensityInBin,'all');
    
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


% --------------- find projections of BinCenters to H, K, L ---------------
MidValue = ([BinCenters; zeros(size(BinCenters)); zeros(size(BinCenters))])'/Matrix + BasicPoint;
outHgrid = MidValue(:,1);
outKgrid = MidValue(:,2);
outLgrid = MidValue(:,3);
% -------------------------------------------------------------------------

% Draw BinCenters points
if Draw == 1
    plot3(outHgrid,outKgrid,outLgrid,'o','markersize',8,'markerfacecolor','b')
end


if ProjectionAxis == 'H'
    BinCoord = outHgrid;
elseif ProjectionAxis == 'K'
    BinCoord = outKgrid;
elseif ProjectionAxis == 'L'
    BinCoord = outLgrid;
elseif ProjectionAxis == 'T'
    BinCoord = BinCenters;
end


end








