function [I, rErr] = rombergDiscrete(intEstimates)
%rombergDiscrete Estimate integrals from a set of initial estimates.
% [I, err]=rombergDiscrete(intEstimates) Estimates an integral using a set of
% precalculated, low accuracy estimates in the column vector
% intEstimates. It returns the estimate of the integral and its relative
% error.

  % Error checking
  if ~isvector(intEstimates)
    error('intEstimates must be a column vector.');
  end

  if ~iscolumn(intEstimates)
    intEstimates = intEstimates';
  end

  intTable = zeros(length(intEstimates));
  relativeErrorTable = zeros(length(intTable));

  intTable(:, 1) = intEstimates;

  for kk = 2:length(intTable)
    for jj = length(intTable)-kk+1:-1:1
      intTable(jj, kk) = ((4^(kk-1))*intTable(jj+1, kk-1) - intTable(jj, kk-1)) / (4^(kk-1) - 1);
    end
  end

  % TODO: Make this work
  for kk = 2:length(intTable)
    relativeErrorTable(kk, 1) = 100 * abs((intTable(1, kk) - intTable(1, kk)) / intTable(1, kk));
  end

  I = intTable(1, end);
  rErr = relativeErrorTable(1, end);
end
