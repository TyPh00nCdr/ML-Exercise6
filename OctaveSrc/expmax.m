## https://de.mathworks.com/help/stats/multivariate-normal-distribution.html
## https://de.mathworks.com/help/stats/mvnrnd.html
## https://de.mathworks.com/help/stats/mvnpdf.html
## https://de.mathworks.com/help/stats/mvncdf.html
## https://octave.sourceforge.io/statistics/function/mvnpdf.html
## https://octave.sourceforge.io/statistics/function/mvncdf.html
## Sigmas: [1.0 0.5; 0.5 1.0], [1, 0.1; 0.1, 0.5], [0.25 0.3; 0.3 1], [0.9 0.4; 0.4 0.3], [1 1.5; 1.5 3]
## [cov (spread) x-axis, rotation; rotation, cov (spread) y-axis] <- must be symmertric i.e. rotation must be equal
## rotation can be negative
## Geometrix Interpretation: https://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
function expmax
  mu1 = [2 2];
  sigma1 = [5 -4; -4 6] / 10;
  
  mu2 = [4 2.5];
  sigma2 = [5 4; 4 6] / 10;
  
  mu3 = [0 1.5];
  sigma3 = [5 4; 4 6] / 10;
  
  cmap = repmat([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], 100*100, 1);
  colormap(cmap);
  # disp(cmap);
  
  [X, Y] = meshgrid (linspace (-1, 5, 100));
  XY = [X(:) Y(:)];
  Z1 = mvnpdf (XY, mu1, sigma1);
  Z1 = reshape (Z1, size (X));
  
  Z2 = mvnpdf (XY, mu2, sigma2);
  Z2 = reshape (Z2, size (X));
  
  Z3 = mvnpdf (XY, mu3, sigma3);
  Z3 = reshape (Z3, size (X));
  contour ([X; X; X], [Y; Y; Y], [Z1; Z2; Z3]);
endfunction
