function parabolicinterp(f::Function, xl::Number, xm::Number, xu::Number, ea::Number, maxIter::Int64)
  # 1D optimisation using parabolic interpolation
  
  xo = NaN; # Current estimate of the optimum
  fxm = f(xm);
  fxu = f(xu);
  fxl = f(xl);

  for i = 1:maxIter
    xoold = xo; # xoxo

    # Solution for xo using simultaneous equations
    fracNumerator = ((xm - xl)^2 * (fxm-fxu)) - ((xm - xu)^2 * (fxm - fxl));
    fracDenominator = ((xm - xl) * (fxm - fxu)) - ((xm - xu) * (fxm - fxl));
    xo = xm - 0.5 * (fracNumerator / fracDenominator);

    if xo < xm
      # The optimum is to the left of xm so make the upper estimate xm and make
      # the middle estimate xo.
      xu = xm;
      fxu = fxm;
      xm = xo;
      fxm = f(xm);
    elseif xo > xm
      # The optimum is to the right of xm so make the lower estimate xm and
      # the middle estimate xo.
      xl = xm;
      fxl = fxm;
      xm = xo;
      fxm = f(xm);
    else
      break;
    end

    # Relative error checks
    if abs((xo - xoold) / xo) * 100 < ea
      break;
    end
  end
  return xo;
end
