function [I] = romberg(f, a, b, nl)
%romberg Numerical Integration of functions using Romberg's method
% I = romberg(f,a,b,nl) Calculates the integral of the function f
% between the limits [a b]. It does so using nl initial estimates of
% the integral using the trapezoid method with successively smaller widths.
% 
% NOTE: f should be a vectorized function such that it can carry out
% element wise multiplication and division on a vector.
%
% NOTE2: nl > 25 should be avoided as they take a _long_ time.

  % Error checks
  if nl <= 0 || mod(nl, 1) ~= 0
    error('Integration level must be a positive integer');
  end

  intTable = nan(nl); % Empty matrix to use as as a divided difference
                      % table
  
  % Calculate the first column of the diff table using the compound trapezium 
  % rule with successively smaller intervals.
  h = (b - a);
  for k = 1:nl
    intTable(k,1) = 0.5 * h *(f(a) + 2*sum(f(a+h:h:b-h)) + f(b));
    h = h / 2;
  end

  % Calculate the remaining elements in the divided difference table using
  % Romberg's method.
  for kk = 2:nl
    for jj = length(intTable)-kk+1:-1:1
      intTable(jj, kk) = ((4^(kk-1) * intTable(jj+1, kk-1)) - intTable(jj, kk-1))/((4^(kk-1)) - 1);
    end
  end
  I = intTable(1,end);
end
