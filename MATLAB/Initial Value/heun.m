function [y] = heun(f, t0, y0, t, h)
%HEUN Solve first order ODEs using Heun's method
%   Using the ODE of the function g, f and the initial conditions t0 and y0, 
%   estimate the value of g(t) using Heun's method with a step size of h.
%   
%   f should be an anonymous function of the form f = @(t,y) ...
%
%
%   NOTE: This may not be a correct implementation of the Heun method. On
%   the last iteration, we use ti+h for the correction step. It's not clear
%   if this should be ti for the last iteration though as ti+h will
%   overshoot.
  y = y0;
  
  if t < t0, h = -abs(h); end
  
  for ti = t0:h:t-h
    yP = y + f(ti, y) * h; % Predictor value
    y = y + 0.5 * h *(f(ti, y) + f(ti+h, yP));
  end
end
