function dist = phasemake(N, M, N_sub)
%%%%% Generate random phase mask

dist = zeros(N);
m = M / N_sub; % segment size
slm_start = (N-M)/2+1; % Internal slm area start point

for i = 1 : N_sub
    for j = 1 : N_sub
        xi = m*(i-1)+slm_start;
        yj = m*(j-1)+slm_start;
        dist(xi : xi+m-1, yj : yj+m-1) = rand*2*pi*8;
    end
end