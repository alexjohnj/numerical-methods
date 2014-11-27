function [I] = simpsonThreeEigthComposite(f, x1, xn, n)
  % simpsonThreeEigthComposite Numerical integration of functions using Simspon's 3/8 composite rule
  %   I = simpsonThreeEigthComposite(f,x1,xn,n) Estimates the integral of f 
  %   between x1 and xn using Simpson's 3/8 rule with n data points. n must
  %   be a multiple of 3.
  if mod(n,3) ~= 0
    error('Must use a multiple of three data points. Please adjust n accordingly');
  end

  h = (xn - x1) / n;
  I = (3/8) * h * (f(x1) + 3*sum(f([x1+h:3*h:xn-2*h])) + ...
  3*sum(f([x1+2*h:3*h:xn-h])) + 2*sum(f([x1+3*h:3*h:xn-3*h])) + f(xn));
end
