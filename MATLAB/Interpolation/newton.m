function [fxint] = newton(x, fx, xint)
  % Interpolate points in a data set
  %
  % fxint = newton(x,fx,xint) Calculates the interpolated value fxint for the 
  % value xint from the vectors x and fx using a Newton interpolating 
  % polynomial.

  % Error checking
  if ~isvector(x) || ~isvector(fx)
    error('X and FX must be column vectors.');
  end

  if ~iscolumn(x)
    x = x';
  end

  if ~iscolumn(fx)
    fx = fx';
  end

  if ~isequal(size(x), size(fx))
    error('X and FX vectors must be the same size.');
  end

  divDifTable = zeros(length(fx)); % Matrix to be used as a divided difference 
                                   % table. 
  divDifTable(:, 1) = fx; % Initialise the first column of the table with fx
  
  for jj = 2:length(divDifTable)
    for ii = length(divDifTable)-jj+1:-1:1
      % Calculate values for each element of the div dif table by 
      % subtracting the elements in the previous column that are in the
      % next row and the same row from each other and then dividing by the
      % difference in the x values.
      divDifTable(ii, jj) = (divDifTable(ii+1, jj-1) - divDifTable(ii, jj-1)) / (x(ii+jj-1, 1) - x(ii, 1));
    end
  end
  
  % Calculate the interpolated value using the div dif table
  fxint = divDifTable(1,1);
  for ii = 2:length(x)
    fxint = fxint + divDifTable(1, ii) * prod(xint - x(1:ii-1, 1));
  end
end
