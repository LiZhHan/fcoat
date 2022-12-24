function [slm1, slm2, df] = slmcoatProA(N_sub, Fs, N, M)
% slm1/slm2 half loading frequency, staggered
% N_sub number of SLM segments
% Fs mark frequency, Fs/2+df : Fs
% N the size of the matrix returned 
% M SLM size

narginchk(4,4);

NN = N_sub*N_sub; % number of segments
df = Fs/NN;
slm1 = zeros(N);slm2 = zeros(N);
m = M / N_sub; % size of segments
slm_start = (N-M)/2+1; % Internal slm area start point

% first half
for i = 1 : N_sub/2
    for j = 1 : N_sub/2
        xi = 2*m*(i-1)+slm_start;
        yj = 2*m*(j-1)+slm_start;
        slm1(xi : xi+m-1, yj : yj+m-1) = NN/2+(i-1)*N_sub+j-1;    % odd row
    end
end

for i = 1 : N_sub/2
    for j = 1 : N_sub/2
        xi = 2*m*(i-1)+slm_start+m;
        yj = 2*m*(j-1)+slm_start+m;
        slm1(xi : xi+m-1, yj : yj+m-1) = NN/2+(i-1)*N_sub+j+N_sub/2-1;    % Even row
    end
end

slm1 = slm1*df;

% second half
for i = 1 : N_sub/2
    for j = 1 : N_sub/2
        xi = 2*m*(i-1)+slm_start; % row
        yj = 2*m*(j-1)+slm_start+m; % column
        slm2(xi : xi+m-1, yj : yj+m-1) = NN/2+(i-1)*N_sub+j-1;    % odd row
    end
end

for i = 1 : N_sub/2
    for j = 1 : N_sub/2
        xi = 2*m*(i-1)+slm_start+m;
        yj = 2*m*(j-1)+slm_start;
        slm2(xi : xi+m-1, yj : yj+m-1) = NN/2+(i-1)*N_sub+j+N_sub/2-1;    % Even row
    end
end
slm2 = slm2*df;