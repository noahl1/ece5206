% filtered back projections file without using imrotate()
clc,clear
ANG=0:0.8:179.2;
load('projs.mat')
PROJECTIONS=p;
PROJECTIONS=PROJECTIONS';
% size of p
GridSize = size(PROJECTIONS,1);


%  FILTER THE PROJECTIONS

[R,C]=size(PROJECTIONS);
%  No Filter applied
% filtPROJECTIONS = PROJECTIONS;

% % make a Ram-Lak filter. it's just abs(f).

% N_ANG = length(ANG);
% N1 = length(R);
% freqs=linspace(-1, 1, N1).';
% myFilter = abs( freqs );
% myFilter = repmat(myFilter, [1 N_ANG]);

% Make Ramp Filter with Cutoff

w = [-pi : (2*pi)/R : pi-(2*pi)/R];
Filt = abs(sin(w));
Filt = Filt(1:513);

% Below we use the Fourier Slice Theorem to filter the image
for i = 1:C
    IMG = fft(PROJECTIONS(:,i));
    FiltIMG = IMG.*Filt';
    %FiltIMG = filter (b, a, IMG);
    FIL(:,i) = ifft(FiltIMG);
end

% Remove any remaining imaginary parts
FIL = real(FIL);
% filtPR = projfilter(PR);

%  Convert the angle to radians and subtract pi so the reconstructed image
%  has the same orientation
 th = pi - ((pi/180)*ANG);

%  Prepare image
 m = length(ANG); 
 BackI = zeros(GridSize,GridSize);

%  Find the middle index of the projections
 midindex = (GridSize+1)/2;

%  Create x and y matrices
 x = 1:GridSize;
 y = 1:GridSize;
 [X,Y] = meshgrid(x,y);
 xprime = X-(GridSize+1)/2;
 yprime = Y-(GridSize+1)/2;

 for i = 1:m
    % Use the backprojection algorithm to determine which areas on the projected
    % images add up 
    filtIndex = round(midindex + xprime*sin(th(i)) - yprime*cos(th(i)));

    % While "in bounds" then add the point
    BackIa = zeros(GridSize,GridSize);
    spota = find((filtIndex > 0) & (filtIndex <= GridSize));
    newfiltIndex = filtIndex(spota);
    BackIa(spota) = FIL(newfiltIndex(:),i);
       
    BackI = BackI + BackIa; 
 end

% BackI = BackI./m;
figure
gg=imshow(BackI,[]);
gg=imadjust(gg,[0 1],[1 0],0.1); %flip values of [0 1] to [1 0], gamma=0 to infinite "contrast"
% imadjust(gg,[0 1],[1 0],0.1)
% figure,imshow(gg)