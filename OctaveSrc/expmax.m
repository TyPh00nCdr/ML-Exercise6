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
  clear X Y XY;
  
  mu1 = [2 2];
  sigma1 = [5 -4; -4 6] / 10;
  cloud1 = mvnrnd (mu1, sigma1, 100);
  
  mu2 = [4 2.5];
  sigma2 = [5 4; 4 6] / 10;
  cloud2 = mvnrnd (mu2, sigma2, 100);
  
  mu3 = [0 1.5];
  sigma3 = [5 4; 4 6] / 10;
  cloud3 = mvnrnd (mu3, sigma3, 100);
  
  scatterplot (cloud1, "r");
  scatterplot (cloud2, "g");
  scatterplot (cloud3, "b");  
  
  contourplot (mu1, sigma1, "r");
  contourplot (mu2, sigma2, "g");
  contourplot (mu3, sigma3, "b");
  
  hold off
  clear X Y XY
endfunction

function phi_j = phi (z, j)
  phi_j =  mean (z == j);
endfunction

function mu_j = mu (z, j, x)
  mu_j = sum ((z == j) * x) / sum ((z == j)); 
endfunction

function mu_j = mu (z, j, x)
  mu_j = sum ((z == j) * x) / sum ((z == j)); 
endfunction

function sigma_j = sigma (z, j, x, mu)
  sigma_j = sum ((z == j) * (x(:) - mu) * (x(:) - mu)', 3) / sum (z == j);  
endfunction

function estimate_weights (z, j, x, mu, sigma, phi)
  
endfunction

function scatterplot (XY, color)
  scatter (XY(:, 1), XY(:, 2), color, "x");
  if (!ishold())   
    hold on
  endif;
endfunction

function contourplot (mu, sigma, color)
  persistent X Y XY
  if isempty(X) || isempty(Y)
    [X, Y] = meshgrid (linspace (-1, 5, 100));
    XY = [X(:) Y(:)];
  endif

  Z = mvnpdf (XY, mu, sigma);
  Z = reshape (Z, size (X));
  contour (X, Y, Z, 1, color);
  if (!ishold())   
    hold on
  endif;
endfunction
