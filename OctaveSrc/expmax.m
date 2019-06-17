## https://de.mathworks.com/help/stats/multivariate-normal-distribution.html
## https://de.mathworks.com/help/stats/mvnrnd.html
## https://de.mathworks.com/help/stats/mvnpdf.html
## https://de.mathworks.com/help/stats/mvncdf.html
## https://octave.sourceforge.io/statistics/function/mvnpdf.html
## https://octave.sourceforge.io/statistics/function/mvncdf.html
## Sigmas: [1.0 0.5; 0.5 1.0], [1, 0.1; 0.1, 0.5], [0.25 0.3; 0.3 1], [0.9 0.4; 0.4 0.3], [1 1.5; 1.5 3]
## [cov (spread) x-axis, rotation; rotation, cov (spread) y-axis] <- must be symmertric i.e. rotation must be equal
## rotation can be negative
function expmax
  mu = [0 0];
  sigma = [.9 -.4; -.4 .3];
  
  [X, Y] = meshgrid (linspace (-3, 3, 100));
  XY = [X(:) Y(:)];
  Z = mvnpdf (XY, mu, sigma);
  Z = reshape (Z, size (X));
  contour (X, Y, Z);
endfunction
