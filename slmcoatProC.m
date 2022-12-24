function [slm_corf,slm_int] = slmcoatProC(N_sub, N, M, I2, slm_cor)

% N the size of the matrix returned
% M SLM size
% N_sub number of SLM segments
% I2 The light intensity sequence obtained at the second half of the desired point
% slm_cor The compensation phase of the previous correction, which is replaced by the result of this time

narginchk(5,5);

slm_int = zeros(N,N);
NN = N_sub*N_sub; % number of segments
slm_corf = slm_cor;
m = M / N_sub; % size of segments
slm_start = (N-M)/2+1; % Internal slm area start point

% The other half replaces the angle calculated in the previous step
F = fft(I2);   %%% Fourier transform
amp = -angle(F);   % angle,(-pi,pi]
amp = amp + (amp<0)*2*pi; % Change the part of (-pi,0) to (pi,2pi)

for i = 1 : N_sub/2
    for j = 1 : N_sub/2
        xi = 2*m*(i-1)+slm_start;
        yj = 2*m*(j-1)+slm_start+m;
        slm_corf(xi : xi+m-1, yj : yj+m-1) = amp(NN/2+(i-1)*N_sub+j);    % odd row
    end
end

for i = 1 : N_sub/2
    for j = 1 : N_sub/2
        xi = 2*m*(i-1)+slm_start+m;
        yj = 2*m*(j-1)+slm_start;
        slm_corf(xi : xi+m-1, yj : yj+m-1) = amp(NN/2+(i-1)*N_sub+j+N_sub/2);    % Even row
    end
end

% binarization
slm_corf = (slm_corf>pi/2&slm_corf<=3*pi/2)*pi;