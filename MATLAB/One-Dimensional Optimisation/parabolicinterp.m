function [xo] = parabolicinterp(f, xl, xm, xu, ea, maxIter)
  % 1D Optimisation using parabolic interpolation.
  %
  %   parabolicinterp(f,xl,xm,xu,maxIter) Estimates the optimum of the
  %   anonymous function f using xl, xm and xu as the lower, middle and
  %   upper guesses for the optimum. It iterates on the estimate until the
  %   relative error falls below ea or maxIter iterations have occurred. 
  
  xo = NaN; % Current estimate of the optimum
  
  % Precompute f(x) for the initial guesses to avoid needlessly
  % recalculating on each iteration
  fxm = f(xm);
  fxu = f(xu);
  fxl = f(xl);
  
  for ii = 1:maxIter
    % Solution for xo using simultaneous equations.
    fracNumerator = ((xm - xl)^2 * (fxm - fxu)) - ((xm-xu)^2*(fxm-fxl));
    fracDenominator = ((xm-xl)*(fxm-fxu)) - ((xm-xu)*(fxm - fxl));
    xlast = xo;
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
      break;
    end
    
    % Calculate relative error and break if it is below our desired
    % accuracy
    if abs((xo - xlast) / xo) * 100 < ea, break; end
  end
end
