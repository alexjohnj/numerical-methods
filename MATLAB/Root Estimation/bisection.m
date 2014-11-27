function [xr] = bisection(f, xl, xu, ea, maxIter)
  %bisection Estimate roots using the bisection method
  %   xr = bisection(f,xl,xu,ea,maxIter) Estimate the root xr of the anonymous
  %   function f inside the search range of xl-xu. Iterate on the estimate
  %   until the relative error of the estimate is less than ea or maxIter 
  %   iterations have taken place.
  
  xr = NaN; % Current estimate of the root
  fxl = f(xl); % Pre-compute fxl to avoid needless function calls
 
  for ii = 1:maxIter
    xrlast = xr;
    xr = 0.5 * (xl+xu);
    fxr = f(xr);
    relativeError = abs((xr - xrlast)/xr) * 100;

    if fxr * fxl == 0
      break;
    elseif fxr * fxl > 0
      xl = xr;
      fxl = fxr;
    else
      xu = xr;
    end
    
    if relativeError < ea
      break;
    end
  end
end
