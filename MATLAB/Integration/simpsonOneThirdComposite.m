function [I] = simpsonOneThirdComposite(f, x1, xn, n)
  % simpsonOneThirdComposite Numerical integration of functions using Simpson's composite 1/3 rule. 
  %
  % I = simpsonOneThirdComposite(f,x1,xn,n) Estimates the integral of the
  % function f between x1 and xn using Simpson's composite 1/3 rule with n 
  % points. n must be odd.
  
  if mod(n,2) == 0
    error('Must use an odd number of points. Please make n odd.');
  end

  h = (xn - x1) / (n - 1);
  I = (1/3) * h * (f(x1) + 4*sum(f([x1+h:2*h:xn-h])) + 2*sum(f([x1+2*h:2*h:xn-2*h])) + f(xn));
end