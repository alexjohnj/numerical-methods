function [xi] = fixedpointiteration(g, xi, ea, maxIter)
%fixedpointiteration Estimate roots using the fixed point iteration method
%   xr = fixedpointiteration(g,xi,ea,maxIter) estimates the root xr of a
%   function using the rearrangement g of the function in terms of x. Using
%   xi as an initial estimate, it iterates on the root until the relative
%   error falls below ea or maxIter iterations occur.   
  for ii = 1:maxIter
    xiold = xi;
    xi = g(xi);
    
    if (abs(xi - xiold) / xi) < ea && ii > 1
      break;
    end
  end
end
