%% Problem 1 - uncertainty principle
% keep track of origin
% in final plot of spectrum the origin is 17th element
% uncertainty principle

%why do i need to keep track of the origin?
%what do i use fftshift for?
%what is proper temporal response plot axis? Microseconds or samples?
%Is the Uncertainty principle 1/(N*dt) or (2*pi)/(N*dt) in this case? How
%to know

clc
clear all
close all

f_n = [-0.7188 0.2242 -0.5723 0.3183 -0.2188 -0.7558 0.1348 -0.6617 0.2812 -1.1617 0.1348 -0.2558 -0.2188 -0.1817 -0.5723 0.7242];
n_n = -8; %starting at n=-8
N = 16;
time_origin = 9;
dt = 0.125;

%Discrete Fourier Transform
F_m = fft(f_n);
%shift
F_m = fftshift(F_m);

%Uncertainty principle
df = 1/(N*dt);
dHz = [-N/2:N/2-1]*df;

tiledlayout(3,1)

nexttile
t = [0 [1:N-1]*dt]-(time_origin*dt);
plot(t,[f_n], '-o')
xlabel('micro seconds') 
ylabel('f[t]') 
title('Temporal Response')

nexttile
plot(dHz,abs(F_m), '-o')
%xlabel('MHz') 
ylabel('Magnitude') 
title('Frequency Response')

nexttile
plot(dHz,angle(F_m), '-o')
xlabel('MHz') 
ylabel('Phase') 


%% Problem 2 - Time Reversal

%For DFT how is reversal theorem adjusted?
%how to explain error to be negligible
%how DFT different from FT
%am i doing this right?

clc
clear all 
close all

pos_n = ([1:2^8]-1);
neg_n = pos_n * -1;


g_n = cos(pos_n) + 1i*sin(pos_n);
plot(pos_n,g_n);

%If g[n] <=> G[m]
G_m = fft(g_n);
recon_g_n = ifft(G_m);
recon_RMSE = sqrt(mean((recon_g_n - g_n).^2));
recon_RMSE
if(recon_RMSE < 0 + 1e-5 + 1i*1e-5)
    fprintf('g[n] <=> G[m] is TRUE');
end

%then y[n]=g[-n] <=> Y[m]=G[-m]
y_n = cos(neg_n) + 1i*sin(neg_n);
Y_m = fft(y_n);
recon_y_n = ifft(Y_m);
recon_RMSE = sqrt(mean((recon_y_n - y_n).^2));
recon_RMSE
if(recon_RMSE < 0 + 1e-5 + 1i*1e-5)
    fprintf('y[n]=g[-n] <=> Y[m]=G[-m] is TRUE');
end


%% Problem 3 - 2d seperability
%if g(x,y)=g_x(x)*g_y(y) then G(wx,wy)=Gx(wx)*Gy(wy),
% where gx(x) <=> Gx(wx) and gy(y) <=> Gy(wy)
clc
clear all
close all
%Can I just prove that the 1d FT of gx(x) and 1d FT of gy(y) can be
%multiplied to get the 2d FT of g(x,y)?
%any function that works?

x = randn(2^8,1);
y = randn(1,2^8);
g_x = x;
g_y = y;
g_x_y = g_x*g_y;

G_wx = fft(g_x);
G_wy = fft(g_y);

Z_wx_wy = G_wx * G_wy;

G_wx_wy = fft2(g_x_y);

recon_RMSE = sqrt(immse(G_wx_wy,Z_wx_wy));
if(recon_RMSE < 0 + 1e-5 + 1i*1e-5)
    recon_RMSE
end
%% Problem 4 - central slice theorem
%CST relates the 1D Fourier transform of the projections to the 2D Fourier transform of the image
clc
clear all
close all
%object domain matrix 2 dimensional: 2^P+1. pad perimeter with 0s
img = double(rgb2gray(imread('butterfly.jpg')))/255;
P = 10; N = 2^P;
img = imresize(img, [N, N]);%Resize
orig_img = img;%to compare
imshow(img); figure
img(end+1, :) = 0; 
img(:, end+1) = 0;
phi = linspace(0,180,1000);
csa = zeros(length(phi), N);
%build matrix of central slices at each angle step, each row is new angle
for p = 1:length(phi)
    %imrotate opposite of angle (clockwise)
    rimg = imrotate(img, p, 'bicubic', 'crop');
    %1d fft for each row -> central slice array
    prj = sum(rimg, 1);
    prj = fftshift(prj(1: N));
    img(:,end) = [];
    prj = fft(prj, N); 
    img(:, end+1) = 0;
    csa(p,:) = fftshift(prj(1:N));
end
plot(csa); figure;
cos_val = cosd(phi);
sin_val = sin(phi);
x = (-N / 2: N / 2-1) .* cosd(phi)';
y = (-N / 2: N / 2-1) .* sind(phi)';
[x_tgt, y_tgt] = meshgrid(-N / 2: N / 2-1, -N / 2: N / 2-1);
img_recon = griddata(x,y,csa,x_tgt,y_tgt);
%Remove NaN
nanloc=find(isnan(img_recon));
img_recon(nanloc) = zeros(size(nanloc));
img_recon=fftshift(img_recon);
img(:,end) = [];
img_recon=ifft2(img_recon);
img(:, end+1) = 0;
img_recon=fftshift(img_recon);
tiledlayout(2,1)
nexttile
imshow(orig_img);
title("Original " + N + "x" + N + " Image");
nexttile
imshow(img_recon);
title("Reconstructed " + N + "x" + N + " Image");

%RMS_error
RMS = sqrt(immse(abs(img_recon),orig_img))





