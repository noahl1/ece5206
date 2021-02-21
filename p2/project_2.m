clear all, clc

projs=load('projs.mat');
projs=projs.p(1:end,1:end-1);

spacing = 1;

%% Apply angle spacing
t=1; %temp index
for i = 1:spacing:size(projs,1)    
    temp(t,:) = projs(i,:);
    t = t + 1;
end
projs=temp;
angle = 0:0.8*spacing:179.2;
angle = angle .*(pi/180);

%% Ram-Lak filter
N = length(projs);
freq = linspace(-1, 1, N);
RamLak = repmat(abs(freq), [size(projs,1) 1 ]);

%% Reconstruction
projs_sh = fftshift(projs,2);
projs_sh_fft = fftshift(fft(projs_sh,[],2),2);
projs_filt = projs_sh_fft .* fftshift(RamLak); %apply filter

% figure
% mesh(abs(projs_filt))

[Wx,Wy] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);
WxKp = -N/2:N/2-1;
Wyk =(sin(-angle)).'*WxKp;
Wxk =(cos(-angle)).'*WxKp;
projs_grid = griddata(Wxk,Wyk,projs_filt,Wx,Wy);

%replacing the nans with zeros for IFFT calc
NanL = find(isnan(projs_grid));
projs_grid(NanL) = zeros(size(NanL));

%figure()
%mesh(abs(projs_grid))

%placing the orgin to the first cell for proper IFFT
projs_circsh = circshift(projs_grid,N/2+1,1);
projs_circsh = circshift(projs_circsh,N/2+1,2);

%shifting back to the center after IFFT and display
projs_ifft = ifft2(projs_circsh);
projs_ifft_sh = fftshift(projs_ifft);

%figure()
%mesh(abs(projs_ifft_sh))

% normalize
projs_norm = abs(projs_ifft_sh) ./ 255;

figure()
imshow(projs_norm)
title("Reconstruction of Image with filter")