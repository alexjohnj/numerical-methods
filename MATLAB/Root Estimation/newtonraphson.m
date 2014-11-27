function [xr] = newtonraphson(f, fprime, xi, ea, maxIter)
%newtonraphson Estimate roots using the Newton-Raphson method
%   xr = newtonraphson(f,fprime,xi,ea,maxIter) estimates the root xr of the
%   anonymous function f using its derivative fprime and an initial
%   estimate of xi. It iterates on the estimate until the relative error
%   falls below ea or maxIter iterations occur.
  xr = xi;
  for ii = 1:maxIter
    xrprev = xr;
    xr = xr - (f(xr)/fprime(xr));
    
    if abs(xr - xrprev)/xr < ea
      break;
    end
  end
end
