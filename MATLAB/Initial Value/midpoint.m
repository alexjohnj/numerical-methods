function [y] = midpoint(f, t0, y0, t, h)
%MIDPOINT Solve first order ODEs using the midpoint method
%   Using the ODE of the function g, f, and the initial conditions t0 and y0, 
%   estimate the value of g(t) using the midpoint method with a step size of h.
%   
%   f should be an anonymous function of the form f = @(t,y) ...
  
  y = y0;
  if t < t0, h = -abs(h); end
  
  for ti = t0:h:t-h
    yHalf = y + 0.5 * h * f(ti, y);
    y = y + f(ti + 0.5*h, yHalf) * h;
  end
end
