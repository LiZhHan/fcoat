function  base = basePro(N, M, ramp,ramk,D)
% According to the given parameters, create a binary partitioned slm basis phase matrix
% N number of matrix elements N*N
% num_area Number of segments num_area*num_area
% ramk Binary Slope Direction

narginchk(5,5);


xx = 0 : M-1; % column
yy = xx; % row
[XX, YY] = meshgrid(xx, yy);
temp = zeros(N);
% In the reverse direction, it should be reduced from the maximum value to 0, so the maximum value is added, and the positive direction is 0
YOFF = (M-1)*(abs(ramk(1))-ramk(1))*ramp/2; 
XOFF = (M-1)*(abs(ramk(2))-ramk(2))*ramp/2;      
temp((N-M)/2+1 : (N-M)/2+M, (N-M)/2+1 : (N-M)/2+M) = (YY*ramk(1) + XX*ramk(2))*ramp + YOFF + XOFF;
base = temp*D/N ;  

               
end