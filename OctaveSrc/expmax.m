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
  clear all;
  
  global ITERATIONS = 50;
  
  mu1 = [2 2];
  sigma1 = [5 -4; -4 6] / 10;
  cloud1 = mvnrnd (mu1, sigma1, 100);
  
  mu2 = [4 2.5];
  sigma2 = [5 4; 4 6] / 10;
  cloud2 = mvnrnd (mu2, sigma2, 100);
  
  mu3 = [0 1.5];
  sigma3 = [5 4; 4 6] / 10;
  cloud3 = mvnrnd (mu3, sigma3, 100);
  
  labels = repelem ([2 3 1], 100).';
  
  figure(1);
  subplot(2, 1, 1);
  scatterplot (cloud1, "r");
  scatterplot (cloud2, "g");
  scatterplot (cloud3, "b");
    
  
  # contourplot (mu1, sigma1, "r");
  # contourplot (mu2, sigma2, "g");
  # contourplot (mu3, sigma3, "b");
  
  global W     = zeros (3, 300);
  global PHI   = [1 / 3, 1 / 3, 1 / 3];
  global MU    = [1 2; 2 2; 3 2];
  global SIGMA = repmat ([1 0; 0 1], [1, 1, 3]);
  X = [cloud1; cloud2; cloud3];
  err = [0];
  
  for i = 1:ITERATIONS
    maximizeparams (X);
    [maxval, c] = max(W, [], 1);
    err(i) = columns(c(c != labels.'));
  end
  
  contourplot (MU(1, :), SIGMA(:, :, 1), "b");
  contourplot (MU(2, :), SIGMA(:, :, 2), "r");
  contourplot (MU(3, :), SIGMA(:, :, 3), "g");
  hold off;
  
  subplot (2, 1, 2);
  plot(1:ITERATIONS, err);
  title("Misclassified Points");
  
  # Reset globals
  MU    = [1 2; 2 2; 3 2];
  SIGMA = repmat ([1 0; 0 1], [1, 1, 3]);
  kmeans (X);
endfunction

function updatephi (x)
  global PHI W
  PHI = mean (W, 2).';
endfunction

function updatemu (x)
  global MU W
  for j = 1:3
    MU(j, :) = (sum(W(j, :)' .* x) ./ sum(W(j, :)))';
  end
endfunction

function updatesigma (x)
  global SIGMA W MU
  for j = 1:3
    y = zeros (size (SIGMA (:, :, j)));
    for i = 1:rows(x)
      y = y + W(j, i) * ((x(i, :) - MU(j, :)).' * (x(i, :) - MU(j, :)));
    end
    SIGMA(:, :, j) = y / sum(W(j, :));
  end
endfunction

function updateweights (x)
  global W MU SIGMA PHI
  s = 0;
  for j = 1:3
    s += mvnpdf(x, MU(j, :), SIGMA(:, :, j)) * PHI(j);
  endfor
  
  for j = 1:3
    W(j, :) = (mvnpdf(x, MU(j, :), SIGMA(:, :, j)) * PHI(j)) ./ s;
  endfor
  c = sum(W);
endfunction

function maximizeparams (x)
  global W MU PHI SIGMA
  updateweights (x);
  updatephi (x);
  updatemu (x)
  updatesigma (x);
endfunction

function kmeans (x)
  global MU ITERATIONS
  c = zeros (300, 1);
  mu = num2cell (MU, 2);
  J = [0];
  labels = repelem ([2 3 1], 100).';
  err = [0];
  
  for loop = 1:ITERATIONS
    for i = 1:rows(x)
      [minVal, c(i)] = min (cellfun (@(mu_j) norm (x(i, :) - mu_j)^2, mu));
    endfor
    for j = 1:size(mu)
      mu(j) = sum ((c(:) == j) .* x) / sum (c(:) == j);
    endfor
    J(loop) = sum(cellfun (@(x_i, c_i) norm (x_i - mu{c_i})^2, num2cell(x, 2), num2cell(c)));
    err(loop) = rows(c(c != labels));
  endfor
  
  figure(2);
  subplot(2, 2, [1 2]);
  scatterplot (x(1:100, :), "r");
  scatterplot (x(101:200, :), "g");
  scatterplot (x(201:300, :), "b"); 
  scatter (cell2mat(mu)(:, 1), cell2mat(mu)(:, 2),  "k", "filled");
  hold off;
  
  subplot(2, 2, 3);
  plot(1:ITERATIONS, J);
  title("Distortion Function");
  
  subplot(2, 2, 4);
  plot(1:ITERATIONS, err);
  title("Misclassified Points");
endfunction

function scatterplot (XY, color)
  scatter (XY(:, 1), XY(:, 2), color, "x");
  if (!ishold())   
    hold on;
  endif;
endfunction

function contourplot (updatemu, sigma, color)
  persistent X Y XY
  if isempty(X) || isempty(Y)
    [X, Y] = meshgrid (linspace (-1, 5, 100));
    XY = [X(:) Y(:)];
  endif

  Z = mvnpdf (XY, updatemu, sigma);
  Z = reshape (Z, size (X));
  contour (X, Y, Z, color);
  if (!ishold())   
    hold on;
  endif;
endfunction
