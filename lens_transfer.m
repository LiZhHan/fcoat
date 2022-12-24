
function U_in = lens_transfer(u_in , k , X, Y, f, Kz, d)
%%% Transfer Function of Light, Simulation of Transmission Through a Lens
% u_in  incident light
% X,Y  horizontal and vertical coordinates of light
% f    lens focal length
% d    Distance to next lens
% k,kz Wave vector and z-axis component


narginchk(7,7);

DL = 125;
u_lenb = u_in .* exp(-1i*k*(X.^2 +Y.^2)/(2*f)).* (X.^2 + Y.^2 <= (DL/2).^2);
U_tempA = fftshift(fft2(u_lenb));
U_tempB = U_tempA .* exp(1i*Kz*d);
U_in = ifft2(ifftshift(U_tempB));