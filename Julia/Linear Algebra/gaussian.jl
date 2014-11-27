function gaussian(A::AbstractMatrix, b::AbstractVector)
  #= Solve linear systems using Gaussian elimination

    Linear systems should be in the form of Ax=b where A is a matrix and b & x 
    are column vectors.
    
    Keyword arguments:
      A -- A matrix containing the coefficients of the variables to solve for
      b -- A column vector containing the RHS of the linear system
    Returns:
      A column vector containing the solutions for x1..xn
  =#
  N, M = size(A);
  x = zeros(M);

  # Copies of A & b to avoid mutation
  cpA = copy(A);
  cpb = copy(b);

  for j = 1:M
    # Maximise the pivot point by looking for larger pivot points in the rows
    # beneath the current point's row.
    maxPivotRowIndex = indmax(abs(cpA[j:end, j])) + (j - 1);
    cpA[maxPivotRowIndex, :], cpA[j, :] = cpA[j, :], cpA[maxPivotRowIndex, :];
    cpb[maxPivotRowIndex, :], cpb[j, :] = cpb[j, :], cpb[maxPivotRowIndex, :];

    # Iterate through the columns then rows of the linear system
    for i = j+1:N
      # If the current element is 0 then it is already reduced
      if cpA[i, j] == 0
        continue
      end

      # Perform row reduction on the current element's row using the pivot
      # point.
      multFactor = cpA[i,j] / cpA[j,j];
      cpA[i, :] = cpA[i, :] - multFactor * cpA[j, :];
      cpb[i] = cpb[i] - multFactor*cpb[j];
    end
  end

  # Perform a back substitution for x1, x2, ..., xn.
  x[end] = cpb[end] / cpA[end, end]; # Calculate xn
  for i = N:-1:1 # Calculate xn-1...x1
    x[i] = (cpb[i] - sum(cpA[i, i+1:N] * x[i+1:N])) / cpA[i,i];
  end
  return x;
end
