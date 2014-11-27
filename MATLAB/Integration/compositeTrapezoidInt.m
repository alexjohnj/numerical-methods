function [I] = compositeTrapezoidInt(f, x1, xn, n)
%compositeTrapezoidInt Integrate functions using the composite trapezoid rule
%   I=compositeTrapezoidInt(f,x1,xn,n) estimates the integral of the
%   function f between x1 through xn using n datapoints with the composite
%   trapezoid rule.
  if n <= 1
    error('Insufficent data points. Try increasing n.');
  end
  
  h = (xn-x1) / (n - 1);
  I = 0.5 * h * (f(x1) + 2*sum(f(x1+h:h:xn-h)) + f(xn));
end
