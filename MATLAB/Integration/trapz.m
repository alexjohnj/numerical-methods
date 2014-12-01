function [I] = trapz(f, a, b, n)
%TRAPZ Integrate functions using the composite trapezoid rule
%   I=trapz(f,a,b,n) estimates the integral of the
%   function f between a through b using n sections with the composite
%   trapezoid rule. Note that f must be vectorized. 
  if n <= 0 || mod(n, 1) ~= 0
    error('n must be a positive integer.');
  end
  
  h = (b-a) / n;
  I = 0.5 * h * (f(a) + 2*sum(f(a+h:h:b-h)) + f(b));
end
