function [] = parabolicinterp(f, xl, xm, xu, ea, maxIter)
  % 1D Optimisation using parabolic interpolation.
  %
  %   parabolicinterp(f,xl,xm,xu,maxIter) Estimates the optimum of the
  %   anonymous function f using xl, xm and xu as the lower, middle and
  %   upper guesses for the optimum. It iterates on the estimate until the
  %   relative error falls below ea or maxIter iterations have occurred. 
  
  % Empty matrix to store successive iterations of the optimum
  results = nan(maxIter, 3);
  
  % Precompute f(x) for the initial guesses to avoid needlessly
  % recalculating on each iteration
  fxm = f(xm);
  fxu = f(xu);
  fxl = f(xl);
  
  for ii = 1:maxIter
    % Solution for xo using simultaneous equations.
    fracNumerator = ((xm - xl)^2 * (fxm - fxu)) - ((xm-xu)^2*(fxm-fxl));
    fracDenominator = ((xm-xl)*(fxm-fxu)) - ((xm-xu)*(fxm - fxl));
    xo = xm - 0.5 * (fracNumerator/fracDenominator); % Current estimate of optimum
    
    if xo < xm
      % The optimum is to the left of xm so make the upper estimate xm, and
      % make the middle estimate xo.
      xu = xm;
      fxu = fxm;
      xm = xo;
      fxm = f(xm);
    elseif xo > xm
      % The optimum is to the right of xm so make the lower estimate xm and
      % the middle estimate xo.
      xl = xm;
      fxl = fxm;
      xm = xo;
      fxm = f(xm);
    else
      results = results(1:ii-1, :);
      break;
    end
    
    % Calculate relative error and break if it is below our desired
    % accuracy
    relativeError = nan(1);
    if ii > 1
      relativeError = abs((xo - results(ii-1, 2)) / xo) * 100;
    end
    
    results(ii, :) = [ii xo relativeError];
    
    if relativeError < ea
      results = results(1:ii, :);
      break;
    end
    
  end
    
  fprintf('Iteration\tOptimum\t\tRelative Error\n');
  fprintf('%d\t\t%.8f\t%.8f\n', results');
end