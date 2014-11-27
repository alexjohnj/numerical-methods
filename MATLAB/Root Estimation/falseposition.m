function [xr] = falseposition(f, xl, xu, ea, maxIter)
  %falseposition Estimate roots using the false position method
  %
  %   falseposition(f,xl,xu,ea,maxIter) Estimates the root of the anonymous
  %   function f between xl-xu. Iterates on the estimate until the relative
  %   error is less than ea or maxIter iterations have occurred.

  xr = NaN; % Current estimate of root
  fxu = f(xu);
  fxl = f(xl);

  for ii = 1:maxIter
    xrold = xr;
    xr = xu - ((fxu*(xl-xu)) / (fxl- fxu));
    fxr = f(xr);
    relativeError = abs((xr - xrold) / xr) * 100;

    if fxr * fxl == 0
      break;
    elseif fxr * fxl > 0
      xl = xr;
      fxl = fxr;
    else
      xu = xr;
      fxu = fxr;
    end

    if relativeError < ea
      break;
    end
  end
end
