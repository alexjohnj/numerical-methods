function [] = golden(f, xl, xu, ea, maxIter)
  % golden One Dimensional Optimisation
  % 
  %   golden(f,xl,xu,ea,maxIter) Estimates the optimum of the anonymous
  %   function f using xl and xu to bracket the search range. It continues
  %   for maxIter iterations or until the approximate error in the estimate
  %   is less than ea.
  k_GOLDEN_RATIO = (1+sqrt(5)) / 2;
  results = zeros(maxIter, 3);
  
  % Calculate x1, x2 & f(x1), f(x2) for the initial lower and upper bounds.
  x1 = xl + (xu - xl)/k_GOLDEN_RATIO; % Upper test point
  x2 = xu - (xu - xl)/k_GOLDEN_RATIO; % Lower test point.
  
  f1 = f(x1);
  f2 = f(x2);
  
  for ii = 1:maxIter
    if f1 > f2
      % The optimum is to the left of x1. Make the upper bound x1, set the
      % upper test point to the previous lower test point and calculate a
      % new lower test point for the new upper bound.
      % 
      xo = x2; % The current estimate of the optimum.
      xu = x1;
      x1 = x2;
      f1 = f2;
      x2 = xu - (xu - xl)/k_GOLDEN_RATIO;
      f2 = f(x2);
      
      err = calculateError(x2, xl, xo);
      results(ii, :) = [ii, xo, err];
      
      if err < ea
        results = results(1:ii, :);
        break;
      end
      
    else
      % The optimum is to the right of x1. Make the lower bound x2, set the
      % lower test point to the previous upper test point and calculate a
      % new upper test point for the new lower bounds.
      xo = x1; % The current estimate of the optimum.
      xl = x2;
      x2 = x1;
      f2 = f1;
      x1 = xl + (xu - xl)/k_GOLDEN_RATIO;
      f1 = f(x1);
      
      err = calculateError(x2, xl, xo);
      results(ii, :) = [ii, xo, err];
      
      if err < ea
        results = results(1:ii, :);
        break;
      end
    end
  end
  fprintf('Iteration\tOptimum\t\tApproximate Error\n');
  fprintf('%d\t\t%.8f\t%.8f\n', results');
  plot(results(:, 2));
end

% Simple function to calculate the approximate error in an estimated 
% optimum value.
function [err] = calculateError(x2, xl, xo)
  err = abs((x2 - xl) / xo) * 100;
end
