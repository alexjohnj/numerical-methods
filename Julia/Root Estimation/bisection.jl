function bisection(f::Function, xl::Number, xu::Number, ea::Number, maxIter::Int64)
  # Estimate roots of functions using the bisection method
  xr = NaN;
  fxl = f(xl);

  for i = 1:maxIter
    xrlast = xr;
    xr = 0.5(xl+xu);
    fxr = f(xr);
    relativeError = abs((xr - xrlast) / xr) * 100;

    if fxr * fxl == 0
      break;
    elseif fxr *fxl > 0
      xl = xr;
      fxl = fxr;
    else
      xu = xr;
    end

    if relativeError < ea
      break;
    end
  end

  return xr;
end
