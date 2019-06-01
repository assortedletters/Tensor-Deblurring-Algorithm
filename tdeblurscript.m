%% initial setup/blurring part
%read image of cat and normalize to be in [0,1]
% A = imread('testcat.jpg');
% M = min(size(A,1),size(A,2));
% A = A(1:M,1:M,:);
% A = double(A);
% A = (A-min(A(:)))/(max(A(:))-min(A(:)));
% A = phantom3d('Modified Shepp-Logan',64);

k = 25; %size of blurring kernel
fro = @(x)sum(sum(sum(x.^2)));



h = zeros(k,k,3);
ctr = size(A,3)/2;

h(:,:,1) = fspecial('gaussian',k,2)/2;
h(:,:,2) = fspecial('gaussian',k,2);
h(:,:,3) = fspecial('gaussian',k,2)/2;
h = h/sum(sum(sum(h)));

hreal = h;

%transforms psf to an optical transfor function (otf) to work in fourier
h = psf2otf(h,size(A));
% h = circshift(h,[0 0 0]);

ht = (double(ttran(h)));
%blurred image
Blurred = ifftn(fftn(A).*h);
% Blurred = circshift(Blurred,[0 0 ctr+1]);

c = 0.01;
E = randn(size(A,1),size(A,2),size(A,3));
E = E/fro(E);
noisyB = Blurred + c*E*fro(Blurred);
imshow(squeeze(Blurred(32,:,:)),[])
Bt = imfilter(A,ht);
%% params
%example params
% Best mse:
% mse = 5.88e-4
% ro = 1e-4;         %TNN regularity param
% gamma_1 = 1e-4;    %Splitting term Param, hor
% gamma_2 = 1e-4;    %Splitting term param, vert
% mu = 1e-4;         %
% tau = 5;          %threshhold paramter

% SUPER sharp edges:
% mse = 6.67e-4
% ro = 1e-6;         %TNN regularity param
% gamma_1 = 1e-4;    %Splitting term Param, hor
% gamma_2 = 1e-4;    %Splitting term param, vert
% mu = 1e-4;         %
% tau = 5;          %threshhold paramter



eps = 1e-2;
tol = 1e-3;
%parameters
ro = 1e-4;         %TNN regularity param
gamma_1 = 1e-1;    %Splitting term Param, hor
gamma_2 = 1e-1;    %Splitting term param, vert
mu = 1e-4;         %
tau = 5;          %threshhold paramter
%initial previous: 0
uprev = zeros(size(Blurred));
%initial estimate
u = noisyB;
T = svt(u,tau);
% T = 0;
%splitting term
[w1 w2] = tgradm(u);
%langragian multiplied: it's a tensor
LAM = -ro*u;
LAMw1 = 0;
LAMw2 = 0;
%Fourier conversion
u = fftn(u); T = fftn(T); LAM = (LAM); B = fftn(noisyB); 
D1 = psf2otf([1 -1], size(u));
D2 = psf2otf([1;-1], size(u));

dif = fro(real(u-uprev));
difprev = inf;
MAX_ITERS = 500;
k = 1;
rel = inf;
myimshow = @(A) imshow(squeeze(A(ctr,:,:)));
while(dif > tol && k < 100 && rel > tol)
    uprev = u;
    %u update:  
    u = (ht.*B + ro*(T-fftn(LAM)) + 0*gamma_1*conj(D1).*fftn(w1-LAMw1/gamma_1) + 0*gamma_2*conj(D2).*fftn(w2-LAMw2/gamma_2)) ./ ... 
        (ht.*h+mu*teye(size(u,1),size(u,3))+0*gamma_1*(D1.*conj(D1))+0*gamma_2*(D2.*conj(D2)) +eps*ones(size(A)));
    %T update
    T = fftn(svt((ifftn(u-fftn(LAM))),tau));
    %w1,w2 update
    gx = real(ifftn(D1.*u)); gy = real(ifftn(D2.*u));
    w1 = sign(gx + LAMw1/gamma_1).*max(abs(gx + LAMw1/gamma_1)-mu/gamma_1,0);
    w2 = sign(gy + LAMw2/gamma_2).*max(abs(gy + LAMw2/gamma_2)-mu/gamma_2,0);
    %Langrangian update
    LAM =  (LAM + ro*(ifftn(u-T)));
    LAMw1 = LAMw1 + gamma_1*(gx-w1);
    LAMw2 = LAMw2 + gamma_2*(gy-w2);
    
    %iteration housekeeping
    k = k+1;
    difprev = dif;
    dif = fro(real(ifftn(double(u-uprev))))
    rel = fro(dif-difprev)/fro(dif);
%     figure
%   imshow(ifftn(w1),[])
%     imshow(ifftn(T))
%     imshow(LAM,[])
    ur = ifftn(u);
%     myimshow(ur)
%     pause(1)
end
u = real(ifftn(u)); T = real(ifftn(T));
% u = circshift(u,[0 0 ctr+1]);
error = immse(u,A);
figure
subplot(1,3,1)
myimshow(u)
title('Restored Image','FontSize',16)
xlabel("mse = "+error)
subplot(1,3,2)
myimshow(noisyB)
title('Degraded Image','FontSize',16)
xlabel("Noise Level = " +100*c +"%")
subplot(1,3,3)
myimshow(A)
title('Ground Truth','FontSize',16)
