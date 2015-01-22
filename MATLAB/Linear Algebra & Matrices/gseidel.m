function [x] = gseidel(A, b, err, maxit)
    % x = gseidel(A,b,err,maxit) produces an estimate of the column vector
    % x where [A]{x}={b} where A is a matrix, b is a column vector, err is
    % the desired minimum error in the estimate and maxit is the maximum
    % iterations on the estimate to perform. The function prints each
    % iteration to the console and stops either when the error in each x
    % estimate falls below err or maxit iterations is reached.

    % Check that matrices have compatible dimensions
    if length(A) ~= length(b)
        error('Incompatible matrix dimensions');
    end

    % Preallocate the x column vector using 0s for the initial estiamtes of
    % xi
    x = zeros(length(b), 1);
    % A vector of 0s that will be used to decide when to stop the
    % iterations. When all the elements = 1 (true), the iterations are
    % stopped.
    stoppingConds = zeros(length(b), 1);
    
    % Loop through iterations
    for iter = 1:maxit
       % Loop through x equations
       for jj = 1:length(A) 
          % Maximise the pivot point by looking for larger pivot points in the
          % rows beneath the current point's row.
          [~, maxPivotRowIndex] = max(abs(A(jj:end, jj)));
          maxPivotRowIndex = maxPivotRowIndex + (jj - 1);
          [A(maxPivotRowIndex, :), A(jj, :)] = deal(A(jj, :), A(maxPivotRowIndex, :));
          [b(maxPivotRowIndex), b(jj)] = deal(b(jj), b(maxPivotRowIndex));
          [x(maxPivotRowIndex), x(jj)] = deal(x(jj), x(maxPivotRowIndex));
          
          xest = b(jj,1); % Temporary variable to store the current estimate of x
          
          % Loop through terms of x equations
          for ii = 1:length(A)
              % Don't include the diagonal in the current row in the
              % estimate
              if ii == jj
                  continue;
              end
              xest = xest - (A(jj,ii) * x(ii,1));
          end
          xest = xest / A(jj, jj);
          errorEstimate = abs(xest - x(jj,1));
          
          if errorEstimate < err
             % If the estimate in the error for the current row is less
             % than our desired accuracy, set the associated row in
             % stoppingConds to true/1. Otherwise set it to false/0
             stoppingConds(jj,1) = true;
             continue; 
          else
              stoppingConds(jj,1) = false;
          end
          
          x(jj, 1) = xest;
       end
       fprintf('%d\t\t%.8f\t%.8f\t%.8f\n', iter, x);
       
       % Check to see if the errors in all the x values have fallen beneath
       % the desired accuracy. Stop if they have.
       if stoppingConds == ones(length(b), 1);
           break;
       end
    end
end
