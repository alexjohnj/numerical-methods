function [I] = richardson(f, a, b, n, n2)
%richardson Numerical integration of functions using Richardson's method
%   I=richardson(f,a,b,n,n2) Integrates the function f through the
%   interval [a b]. It calculates two initial estimates of the integral
%   using the compound trapezium rule with n and n2 intervals
%   respectively. It then uses Richardson's method to produce an improved
%   estimate of the integral.

  if n <= 0 || n2 <= 0 || mod(n,1) ~= 0 || mod(n2, 1) ~= 0
    error('The number of intervals must be a positive integer');
  end

  [I1, h1] = ctrap(f, a, b, n);
  [I2, h2] = ctrap(f, a, b, n2);
  
  I = I2 + ((I1 - I2) / (1 - (h1/h2)^2));
end
