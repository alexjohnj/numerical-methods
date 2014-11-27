function [I] = trapezoidInt(f, x1, x2)
%trapezoidInt Integrate functions using the trapezoid rule
%   I=trapezoidInt(f,x1,x2) estimates the integral I of the function f
%   between the limits x1 and x2 using a single trapezoid.
  I = (x2 - x1) * 0.5 * (f(x1) + f(x2));
end
