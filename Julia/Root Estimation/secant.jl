function secant(f::Function, xi::Number, dx::Number, ea::Number, maxIterations::Integer)
  # Estimate the root of functions using the secant method. 
  # 
  # Keyword arguments:
  #  f -- Anonymous function to find the root of
  #  xi -- Initial estimate for the root
  #  dx -- Pertubation to use for estimating df/dx
  #  ea -- Minimum relative error for a solution
  #  maxIterations -- Maximum number of iterations to carry out
  #
  # Returns:
  #  An estimate of the root of the function

  xo = xi; # Current estimate of root

  for i = 1:maxIterations
    xoLast = xo;
    xo = xo - (f(xo) * dx) / (f(xo+dx) - f(xo));
    relativeError = abs((xo - xoLast) / xo) * 100;

    (relativeError < ea || isnan(relativeError)) && break;
  end

  return xo;
end
