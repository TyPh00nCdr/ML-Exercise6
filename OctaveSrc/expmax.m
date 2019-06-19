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
  cloud1 = mvnrnd (mu1, sigma1, 100);
  
  mu2 = [4 2.5];
  sigma2 = [5 4; 4 6] / 10;
  cloud2 = mvnrnd (mu2, sigma2, 100);
  
  mu3 = [0 1.5];
  sigma3 = [5 4; 4 6] / 10;
  cloud3 = mvnrnd (mu3, sigma3, 100);
  
  ax = axes ();
  scatter (ax, cloud1(:, 1), cloud1(:, 2), "r", "x");
  
  hold on;
  
  scatter (ax, cloud2(:, 1), cloud2(:, 2), "g", "x");
  scatter (ax, cloud3(:, 1), cloud3(:, 2), "b", "x");
  
  
  
  [X, Y] = meshgrid (linspace (-1, 5, 100));
  XY = [X(:) Y(:)];
  Z1 = mvnpdf (XY, mu1, sigma1);
  Z1 = reshape (Z1, size (X));
  contour (X, Y, Z1, 1, "r");
  
  Z2 = mvnpdf (XY, mu2, sigma2);
  Z2 = reshape (Z2, size (X));
  contour (X, Y, Z2, 1, "g");
  
  Z3 = mvnpdf (XY, mu3, sigma3);
  Z3 = reshape (Z3, size (X));
  contour (X, Y, Z3, 1, "b");
  
  hold off;
endfunction
