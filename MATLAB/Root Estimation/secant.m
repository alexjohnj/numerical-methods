function [xo] = secant(f, xi, dx, ea, maxIt)
% secant Estimates roots of functions using the secant method
%   xo=secant(f,xi,dx,ea,maxIt) Estimates the root of the function f using
%   the intial estimate xi and a pertubation dx. Iterates on the estimate
%   until the relative error between estimates falls below ea or maxIt
%   iterations have been carried out.
  xo = xi;
  for ii = 1:maxIt
    xoLast = xo;
    xo = xo - (f(xo) * dx) / (f(xo+dx) - f(xo));

    relativeError = abs((xo - xoLast) / xo) * 100;

    if relativeError < ea || isnan(relativeError)
      break;
    end
  end
end
