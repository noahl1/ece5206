clear;
clc;
P=load('projs.mat');
pnew=P.p(1:end,1:end-1);
angle=0:0.8*(pi/180):179.2*(pi/180);
a=1;
for i=1:5:size(pnew,1)
    
pnewNew(a,:)=pnew(i,:);
a=a+1;
end
pnew=pnewNew;
angle=0:0.8*(pi/180)*5:179.2*(pi/180);
%% Ram-Lak filter
N=length(pnew);
FreqLength=linspace(-1, 1, N).';
RamLak=abs( FreqLength' );

% NanL=find(RamLak<0.1);
% RamLak(NanL)=ones(size(NanL))*.1;


RamLak=repmat(RamLak, [size(pnew,1) 1 ]);

%% reconstructing 

SinoNewshift=fftshift(pnew,2);
SinoNewshiftFFT=fftshift(fft(SinoNewshift,[],2),2).*fftshift(RamLak);%putting the filter in here
figure()
mesh(abs(SinoNewshiftFFT))
[Wx,Wy]=meshgrid(-length(SinoNewshift)/2:length(SinoNewshift)/2-1, -length(SinoNewshift)/2:length(SinoNewshift)/2-1);
WxKp=-length(SinoNewshift)/2:length(SinoNewshift)/2-1;
Wyk=(sin(-angle)).'*WxKp;
Wxk=(cos(-angle)).'*WxKp;
SinoNewshiftFFTGrid=griddata(Wxk,Wyk,SinoNewshiftFFT,Wx,Wy);
%replacing the nans with zeros for IFFT calc
NanL=find(isnan(SinoNewshiftFFTGrid));
SinoNewshiftFFTGrid(NanL)=zeros(size(NanL));
figure()
mesh(abs(SinoNewshiftFFTGrid))
%placing the orgin to the first cell for proper IFFT
SinoNewshiftFFTGrid=circshift(SinoNewshiftFFTGrid,length(SinoNewshiftFFTGrid)/2+1,1);
SinoNewshiftFFTGrid=circshift(SinoNewshiftFFTGrid,length(SinoNewshiftFFTGrid)/2+1,2);
%shifting back to the center after IFFT and display
SinoIFFT=ifft2(SinoNewshiftFFTGrid);
SinoIFFT=fftshift(SinoIFFT);
figure()
mesh(abs(SinoIFFT))
% normalizing so I can use imshow and so i can correlate
Max = max(max(abs(SinoIFFT)));
SinoIFFTNrom=abs(SinoIFFT)./abs(Max);

figure()
imshow(SinoIFFTNrom)
title("Reconstruction of Image with filter")