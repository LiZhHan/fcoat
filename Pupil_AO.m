%% initial parameter settings
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   initial parameter settings   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamda = 589e-3; % wavelength
F = 625;  %Focal length
D = 500;  %Lens diameter

M = 256; %DMD or SLM pixels
s = 4;  %FFT expansion multiple
N = s * M; %total matrix,expand prevent FFT superposition
DL = D / s ; %Beam diameter, s times smaller than the lens
N_slmsub = 32; %SLM segments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   coordinate settings   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  object coordinates %%%
dx1 = D / N;
x1 = -D/2 : dx1: D/2-dx1 ; y1 = x1;
[X1,Y1] = meshgrid(x1, x1);
%%% Interpolation parameter settings %%%
%%%In order to draw the light intensity distribution, find the full width at half maximum
Num = 50;
sx_interp =  -dx1*Num:dx1/20:dx1*Num-dx1/30;% Interpolation number
[SX_interp,SY_interp] = meshgrid(sx_interp,sx_interp);%Interpolation frequency
%%%  Image coordinates  %%%
%%%  Pixel coordinates  %%%
label = DL/2 + 20*dx1;
label1 = DL - 100*dx1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  wave vector setting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 2*pi/lamda;
kx = linspace(-pi/dx1, pi/dx1, length(x1)+1);
kx = kx(1 : end-1);
[Kx,Ky] = meshgrid(kx, kx);
Kz = sqrt(k^2-Kx.^2-Ky.^2); %wave vector sampling; spectrum sampling
clear kx Kx Ky
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_in = ones(N, N) .* (X1.^2 + Y1.^2 <= (DL/2).^2);%light source

%% Ideal imaging, through a lens
U_temp1 = fftshift(fft2(u_in));
U_temp2 = U_temp1.* exp(1i*k*F);
u_in1 = ifft2(ifftshift(U_temp2));
u_im = lens_transfer(u_in1 , k , X1, Y1, F, Kz, F);

IS_id = abs(u_im) .^2;
IS_ids = IS_id / sum(sum(IS_id));
IS_mid = IS_ids(N/2+1,N/2+1);
select_data_interp = interp2(X1,Y1,IS_ids,SX_interp,SY_interp,'bicubic');

figure(1);imagesc(x1, y1, IS_ids);axis square;colormap(hot);axis([-label, label, -label, label]);
colorbar; title('Ideal Airy','fontsize',14);xlabel('X/\mum','fontsize',14); ylabel('Y/\mum','fontsize',14);set(gca,'fontSize',14);
% After interpolation
figure(2);imagesc(sx_interp,sx_interp, select_data_interp);axis square;colormap(hot);
colorbar; title(['Ideal Airy : I=',num2str(IS_mid),'(a.u.)'],'fontsize',14);xlabel('X/\mum','fontsize',14); ylabel('Y/\mum','fontsize',14);set(gca,'fontSize',14);

% Gaussian fitting to find full width at half maximum
figure(3);
I_id_midy = select_data_interp(:,20*Num+1);
fitt = fit(sx_interp.', I_id_midy, 'gauss1');
FMWHX = 2*sqrt(log(2))*fitt.c1;
subplot(1,2,1);plot(fitt,sx_interp,I_id_midy);axis([-label, label,0,0.02]);title(['X-axis : FMWH=',num2str(FMWHX),'\mum'],'fontsize',14);
xlabel('x/\mum','fontsize',14); ylabel('Intensity(a.u.)','fontsize',14);set(gca,'fontSize',14);

I_id_midx = select_data_interp(20*Num+1,:);
fitt = fit(sx_interp.', I_id_midx.', 'gauss1');
FMWHY = 2*sqrt(log(2))*fitt.c1;
subplot(1,2,2);plot(fitt,sx_interp,I_id_midy);axis([-label, label,0,0.02]);title(['Y-axis : FMWH=',num2str(FMWHY),'\mum'],'fontsize',14);
xlabel('Y/\mum','fontsize',14); ylabel('Intensity(a.u.)','fontsize',14);set(gca,'fontSize',14);
%% scattering
phase_mask = phasemake(N, N, N_slmsub*4);% random generation of scattering medium

u_i = lens_transfer(u_in1 ,k , X1, Y1, F, Kz, F/2);
u_temp = u_i .* exp(1i*phase_mask);
U_tempA = fftshift(fft2(u_temp));
U_tempB = U_tempA .* exp(1i*Kz*F/2);
u_im = ifft2(ifftshift(U_tempB));
I_s = abs(u_im) .^2;
I_s = I_s / sum(sum(I_s));
I_S_mid = I_s(N/2+1,N/2+1);
Scattering_index = -log(I_S_mid/IS_mid);% Find the scattering coefficient
figure(4);imagesc(x1, y1, phase_mask);axis square;colormap(hot);axis square; title('Phase-mask','fontsize',14);
colorbar;xlabel('X/\mum','fontsize',14); ylabel('Y/\mum','fontsize',14); set(gca,'fontSize',14);

I_s_interp = interp2(X1,Y1,I_s,SX_interp,SY_interp,'bicubic');

figure(5);imagesc(sx_interp,sx_interp, I_s_interp);axis square;colormap(hot);title(['Scattering : I=',num2str(I_S_mid),'(a.u.)'],'fontsize',14);
colorbar;xlabel('X/\mum','fontsize',14); ylabel('Y/\mum','fontsize',14); set(gca,'fontSize',14);

%%  Pupil-AO Compensation

tic

loopnum = 1;% number of iterations
NFFT = N_slmsub*N_slmsub*2;
Fs = 6.4; % The marker frequency
[slm_coatA, slm_coatB, df] = slmcoatProA(N_slmsub, Fs, N, M); % Returns the staggered frequency distribution
I1 = zeros(1, NFFT); % Intensity of the first half desired point
I2 = zeros(1, NFFT); % Intensity of the second half desired point
slm_corn = zeros(N); % The compensation phase obtained from the previous calibration, initially set to 0
itnum = 0;

for loop = 1:loopnum % number of iterations
    
    slm_corf = zeros(N); % The compensation phase obtained by the first half of the correction, set to 0 every cycle
    num_I = 0; % count
    
    % first half
    for t = (0 : NFFT-1) / (Fs*2)   %%%% time
        
        t1=t*2*Fs;
        % binarization
        slm_BcoatA=(mod(slm_coatA*t,1)>=1/4&mod(slm_coatA*t,1)<3/4)*0.5;
        
        u_SLM2 = u_in.* exp(1i*slm_corn) .* exp(1i*2*pi*slm_BcoatA);    % PAO modulation
        U_SLM2temp1 = fftshift(fft2(u_SLM2));
        U_SLM2temp2 = U_SLM2temp1.* exp(1i*k*F);
        u_in1 = ifft2(ifftshift(U_SLM2temp2));
        u_i = lens_transfer(u_in1 , k , X1, Y1, F, Kz, F/2);
        u_temp = u_i .* exp(1i*phase_mask);
        U_tempA = fftshift(fft2(u_temp));
        U_tempB = U_tempA .* exp(1i*Kz*F/2);
        u_focus1 = ifft2(ifftshift(U_tempB));
        
        
        ISC1 = abs(u_focus1) .^2;
        num_I = num_I +1;
        I1(num_I) = ISC1(N/2+1, N/2+1);% Feedback signal
        itnum = itnum+1;
        I3(itnum) = ISC1(N/2+1, N/2+1);
        
    end
    
    % assign half
    slm_corf = slmcoatProB(N_slmsub, N, M, I1, slm_corf); % Update the compensation phase and replace the result of this time
    num_I = 0; % count
    phase_1 = mod(slm_corf, 2*pi);
    
    % the other half
    for t = (0 : NFFT-1) / (Fs*2)
        
        t1=t*2*Fs;
        % binarization
        slm_BcoatB=(mod(slm_coatB*t,1)>=1/4&mod(slm_coatB*t,1)<3/4)*0.5;
        
        u_SLM2 = u_in.* exp(1i*(slm_corn+slm_corf)) .* exp(1i*2*pi*slm_BcoatB);
        U_SLM2temp1 = fftshift(fft2(u_SLM2));
        U_SLM2temp2 = U_SLM2temp1.* exp(1i*k*F);
        u_in1 = ifft2(ifftshift(U_SLM2temp2));
        u_i = lens_transfer(u_in1 , k , X1, Y1, F, Kz, F/2);
        u_temp = u_i .* exp(1i*phase_mask);
        U_tempA = fftshift(fft2(u_temp));
        U_tempB = U_tempA .* exp(1i*Kz*F/2);
        u_focus2 = ifft2(ifftshift(U_tempB));
        
        ISC2 = abs(u_focus2) .^2;
        num_I = num_I +1;
        I2(num_I) = ISC2(N/2+1, N/2+1);
        itnum = itnum+1;
        I3(itnum) = ISC2(N/2+1, N/2+1);
        
    end
    slm_corf = slmcoatProC(N_slmsub, N, M, I2, slm_corf);
    slm_corn = mod(slm_corn+slm_corf, 2*pi);
end
slm_corn1 = slm_corn;
phase_2 = slm_corn1;
toc
%% compensation results after nth iteration

u_SLM2 = u_in .* exp(1i*slm_corn1);
U_SLM2temp1 = fftshift(fft2(u_SLM2));
U_SLM2temp2 = U_SLM2temp1.* exp(1i*k*F);
u_in1 = ifft2(ifftshift(U_SLM2temp2));
u_i = lens_transfer(u_in1 , k , X1, Y1, F, Kz, F/2);
u_temp = u_i .* exp(1i*phase_mask);
U_tempA = fftshift(fft2(u_temp));
U_tempB = U_tempA .* exp(1i*Kz*F/2);
u_focus = ifft2(ifftshift(U_tempB));

I_cor = abs(u_focus) .^2;
itnum = itnum+1;
I3(itnum) = I_cor(N/2+1, N/2+1);
I_cor = I_cor / sum(sum(I_cor));
I_cor_interp = interp2(X1,Y1,I_cor,SX_interp,SY_interp,'bicubic');
I_cor_mid = I_cor(N/2+1,N/2+1);

figure(6);imagesc(sx_interp,sx_interp, I_cor_interp);axis square;colormap(hot);title(['PAO : I=',num2str(I_cor_mid),'(a.u.)'],'fontsize',14);
colorbar;xlabel('X/\mum','fontsize',14); ylabel('Y/\mum','fontsize',14); set(gca,'fontSize',14);

figure(7);
I_id_midy =I_cor_interp(:,20*Num+1);
fitt = fit(sx_interp.', I_id_midy, 'gauss1');
FMWHX = 2*sqrt(log(2))*fitt.c1;
subplot(1,2,1);plot(fitt,sx_interp,I_id_midy);axis([-label, label,0,0.013]);title(['X-axis : FMWH=',num2str(FMWHX),'\mum'],'fontsize',14);
xlabel('x/\mum','fontsize',14); ylabel('Intensity(a.u.)','fontsize',14);set(gca,'fontSize',14);

I_id_midx = I_cor_interp(20*Num+1,:);
fitt = fit(sx_interp.', I_id_midx.', 'gauss1');
FMWHY = 2*sqrt(log(2))*fitt.c1;
subplot(1,2,2);plot(fitt,sx_interp,I_id_midy);axis([-label, label,0,0.013]);title(['Y-axis : FMWH=',num2str(FMWHY),'\mum'],'fontsize',14);
xlabel('Y/\mum','fontsize',14); ylabel('Intensity(a.u.)','fontsize',14);set(gca,'fontSize',14);