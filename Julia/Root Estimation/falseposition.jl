function falseposition(f::Function, xl::Number, xu::Number, ea::Number, maxIter::Int64)
  # Estimate roots of functions using the false position method
  xr = NaN; # Initial estimate of the root
  fxu = f(xu);
  fxl = f(xl);

  for i = 1:maxIter 
    xrold = xr;
    xr = xu - fxu*(xl-xu) / (fxl-fxu);
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
  return xr;
end
