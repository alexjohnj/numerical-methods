function [y] = euler(f, t0, y0, t, h)
%EULER Solve first order ODEs using Euler's method
%   Using the ODE of the function g, f and the initial conditions t0 and y0, 
%   estimate the value of g(t) using Euler's method with a step size of h.
%   
%   f should be an anonymous function of the form f = @(t,y) ...


  if t < t0, h = -abs(h); end

  y = y0;
  for ti = t0:h:t-h
    y = y + f(ti, y)*h; 
  end
end
