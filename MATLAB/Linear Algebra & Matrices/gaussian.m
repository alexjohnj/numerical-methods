function [X] = gaussian(A, b)
  % Solve a linear system using Gaussian elimination
  %
  % X = gaussian(A,b) Given that A is a square matrix and b is a column vector,
  % find the column vector X that satisfies the equation AX=b.
  [N, M] = size(A);
  X = zeros(M, 1);
  
  for jj = 1:M
    % Maximise the pivot point by looking for larger pivot points in the
    % rows beneath the current point's row.
    [~, maxPivotRowIndex] = max(abs(A(jj:end, jj)));
    maxPivotRowIndex = maxPivotRowIndex + (jj - 1);
    [A(maxPivotRowIndex, :), A(jj, :)] = deal(A(jj, :), A(maxPivotRowIndex, :));
    [b(maxPivotRowIndex), b(jj)] = deal(b(jj), b(maxPivotRowIndex));
    
    % Iterate through the columns of the linear system, then each row in
    % said column.
    for ii = jj+1:N
      % If the current element is already 0, it is reduced, so move on to
      % the next one.
      if A(ii, jj) == 0
        continue;        
      end
      
      % Perform row reduction on the current element's row using
      % the pivot point.
      multFactor = (A(ii, jj) / A(jj, jj));
      A(ii, :) = A(ii, :) - (multFactor * A(jj, :));
      b(ii) = b(ii) - multFactor*b(jj);
    end
  end
  % Finally, perform a back subsitution to solve for x1, x2, ..., xn
  X(end) = b(end) / A(end, end); % Calculate xn
  
  % Calculate xn-1 ... x1
  for ii = N-1:-1:1
    summation = 0;
    for jj = ii+1:N
      summation = summation + (A(ii, jj) * X(jj));
    end
    X(ii) = (b(ii) - summation) / A(ii, ii);
  end
end
