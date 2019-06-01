%% setup
% 
% num = 6;
% image_name = cell(1,num);
% image_names = cell(1,num);
% ReErr = cell(1,num);
% PSNRv = cell(1,num);
% image_name{1} ='2.png';
% image_name{2} ='circle.png';
% image_name{3} ='toy.png';
% image_name{4} ='car4.png';
% image_name{5} ='lena256.png';
% image_name{6} ='letters6.png';
% 
% image_names{1} ='2s.png';
% image_names{2} ='circles.png';
% image_names{3} ='toys.png';
% image_names{4} ='car4s.png';
% image_names{5} ='lena256s.png';
% image_names{6} ='letters6s.png';
% 
% I1 = double(imread(image_name{ii}))./255;
% I = I1(:,:,1);
% A = I;
A = phantom(256);
h = fspecial('gaussian', 5,2);
B = imfilter(A,h,'circular','conv');
mask = zeros(size(B));
mid = length(mask)/2;
step = mid/2;
mask(mid-step:mid+step,mid-step:mid+step) = 1;
mask = ones(size(B));
Is = activecontour(B,mask);
Is = sign(B);
Is = mask;
% imshow(B)








%% Blind Deblur

%initialization
u = zeros(5); u(1) = 1; u = fftshift(u);
[mm, nn] = size(u);
s = mm*nn;
g = B;
fb = 0;

v = B*0;
wx = 0; wy = 0;
vbar = 0; ubar = u; wxbar = 0; wybar = 0;

px = Dx(B)*0; py = Dy(B)*0;
q = 0; r = u;

Tu = 4e-2;
Tv = 6e-1;
Tw = 4e-2;

Tp = 2e-2; 
Tq = 2e-2; 
Tr = 2e-2;

alpha = 1e3;
beta = 1e-3;
gamma = 3e-3;
theta= 1;

maxIters = 100;
k = 0;
uprev = u;
    vprev = v;
    wxprev = wx;
    wyprev = wy;
    
while(k < maxIters) 
    k = k+1;
    %p update
    px = px+gamma*Tp*Dx(vbar);
    py = py+gamma*Tp*Dy(vbar);
    normp = sqrt(px.*px+py.*py);
    px = px./max(1,normp);
    py = py./max(1,normp);

    % q update
    Gubar = imfilter(B,ubar,'circular','same','conv');
    q = q+Tq*(vbar - Gubar);
    
    %r update
    r = r + Tr*(ubar+ Dxt(wxbar) + Dyt(wybar));
    
    %u update
    GTq = real(otf2psf(conj(fft2(B)).*fft2(q),size(u)));
    rhs = uprev + Tu*alpha+Tu*(GTq-r);
    u = rhs - Tu*alpha/(Tu*alpha*s+1)*sum(rhs(:));
    u_sym=0.5*(u+u');   %symmetric
    u=0.5*(u_sym+rot90(u_sym,2)');  %symmetric and asymmetric
    
    
    %v update
        tmp_mask=(1-Is)+Is.*double(v<0);
        vup=vprev+Tv*(1-Is).*B-Tv*(gamma*Dxt(px)+gamma*Dyt(py)+q);
        vdown=1+Tv*tmp_mask;
        v=vup./vdown;
    
    
    % w update
    powx = wx+Tw*(Dx(r));
    powy = wy + Tw*(Dy(r));
    [tempx, tempy] = Project(powx,powy,beta*Tw,1);
    wx = powx-tempx;
    wy = powy-tempy;
    
    % bars
    % update ubar,vbar
    ubar = u+theta*(u-uprev);
    vbar = v+theta*(v-vprev);
    wxbar= wx+theta*(wx-wxprev);
    wybar= wy+theta*(wy-wyprev);
    
    % prev updates
    uprev = u;
    vprev = v;
    wxprev = wx;
    wyprev = wy;
    
    
    
    deblured_img=v.*Is;
    deblured_img=min(max(0,deblured_img),1);
    deblured_img = (1-Is).*A+Is.*deblured_img;

end
X = imfilter(B,u,'circular','same','conv');
imshow(X);
figure; imshow(deblured_img)


%% y only u work?

  u = fspecial('gaussian',7,1e-10);% inital inverse psf = delta
    [mm,nn] = size(u);L = mm*nn; e = ones(L,1);
    v = B; q = 0*I; r = 0*u;
    px = 0*I; py = 0*I;wx = 0*u; wy= 0*u;
    u_bar = u; v_bar = v; u_old = u; v_old = v;
    wx_bar=wx; wy_bar=wy; wx_old=wx; wy_old=wy;
    t_u = 0.04;t_v = 0.6; t_w=0.04;
    t_p = 0.02; t_q = 0.02; t_r=0.02;
    theta = 1; gamma=0.003; alpha = 1000; beta=0.001;
    k = 0; crit = 1;
    t0 = cputime;
    
    % while crit>1e-4
    for i=1:1000
        k=k+1;
        
        % update p
        px = px+gamma*t_p*Dx(v_bar);
        py = py+gamma*t_p*Dy(v_bar);
        normp = sqrt(px.*px+py.*py);
        px = px./max(1,normp);
        py = py./max(1,normp);
        
        % update q
        Gubar = imfilter(B,u_bar,'circular','same','conv');
        q = q+t_q*(v_bar-Gubar);
        
        % update r
        r=r+t_r*(u_bar+Dxt(wx_bar)+Dyt(wy_bar));
        
        % update u
        GTq = real(otf2psf(conj(fft2(B)).*fft2(q),size(u)));
        rhs = u_old+t_u*alpha+t_u*(GTq-r);
        u = rhs-t_u*alpha/(t_u*alpha*L+1)*sum(rhs(:));
        
        u_sym=0.5*(u+u');   %symmetric
        u=0.5*(u_sym+rot90(u_sym,2)');  %symmetric and asymmetric
%             u = u./sum(u(:));
        
        % update v
        tmp_mask=(1-mask_S)+mask_S.*double(v<0);
        vup=v_old+t_v*(1-mask_S).*I-t_v*(gamma*Dxt(px)+gamma*Dyt(py)+q);
        vdown=1+t_v*tmp_mask;
        v=vup./vdown;
        %     v=min(max(0,v),1);
        
        % update w
        tepx= wx_old-t_w*Dx(r);
        tepy= wy_old-t_w*Dy(r);
        [powx,powy] = Project(tepx,tepy,beta*t_w,2);
        wx = tepx - powx;
        wy = tepy - powy;
        
        
        % update ubar,vbar
        u_bar = u+theta*(u-u_old);
        v_bar = v+theta*(v-v_old);
        wx_bar= wx+theta*(wx-wx_old);
        wy_bar= wy+theta*(wy-wy_old);
        
        
        Bu = imfilter(B,u,'circular','same','conv');
        
        %     deblured_img=deblured_img.*mask_S/max(max(deblured_img.*mask_S));  % scaling
        deblured_img=v.*mask_S;
        deblured_img=min(max(0,deblured_img),1);
        
        sse = sum(sum((I*255-deblured_img*255).^2.*mask_S));
        mse = sse/(sum(sum(mask_S)));
        mypsnr(k) = 10*log10(((2^8-1)^2)/mse); % for support image
        
        crit = norm(v-v_old,'fro')./norm(v,'fro');
        
        critt(k)=crit;
        if crit < 1*1e-5;
%             break;
        end
        u_old = u;
        v_old = v;
        wx_old = wx;
        wy_old = wy;
        t1 = cputime - t0;
        deblured_img = (1-mask_S).*B+mask_S.*deblured_img;
        %     imagesc(deblured_img);colormap gray;
        %     title(sprintf('PSNR:%5.2fdB. Time:%3.1fs. Iter:%d',mypsnr(k), t1, k));
        %     drawnow
    end