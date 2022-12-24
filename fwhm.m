function width = fwhm(x,y,N)

% function width = fwhm(x,y)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% Rev 1.2, April 2006 (Patrick Egan)

y_max = max(max(y));
nx = find(y==y_max);
start = interp1(y(1 : nx),x(1 : nx),0.5*y_max,'pchip') ;
width = interp1(y(nx+1 : N),x(nx+1 : N),0.5*y_max,'pchip') - start ;
width = width ;% unit um

