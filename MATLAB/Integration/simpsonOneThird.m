function [I] = simpsonOneThird(f, x1, x2, x3)
  % simpsonOneThird Numerical integration of functions using Simpson's 1/3 rule.
  %   I = simpsonOneThird(f,x1,x2,x3) estimates the integral of f between x1 
  %   and x3 uing Simpson's 1/3 rule with a midpoint of x2.

  I = (x3 - x1) * (1/6) * (f(x1) + 4*f(x2) + f(x3));
end
