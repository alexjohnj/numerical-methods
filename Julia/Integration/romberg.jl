function romberg(f::Function, a::Number, b::Number, nl::Int)
  #= Numerical Integration of functions using Romberg's method.
  Keyword arguments:
  f -- Anonymous functon to integrate
  a -- Lower limit of integration
  b -- Upper limit of integration
  nl -- The level of integration to go to.

  Returns:
  The result of the integration of f between a & b
  =#

  intTable = zeros(nl, nl); # Empty matrix to use as div diff table

  # Calculate the first column of the div diff table using the trapezium rule.
  h = (b - a);
  for k = 1:nl
    intTable[k,1] = 0.5 * h * (f(a) + 2*sum(f([a+h:h:b-h])) + f(b));
    h /= 2;
  end

  # Calculate the remaining elements using Romberg's method.
  for k = 2:nl, j = nl-k+1:-1:1
    intTable[j, k] = ((4^(k-1))intTable[j+1, k-1] - intTable[j, k-1])/((4^(k-1)) - 1);    
  end

  return intTable[1, end];
end

function romberg(intEstimates::Array{Float64,1})
  #= Estiamte integrals from sets of data using the Romberg Method.

  Keyword arguments:
  intEstimates -- An array of initial estimates of integrals

  Returns:
    A tuple containing the integral and the relative error in the estimate 
  =# 
  intTable = zeros(length(intEstimates), length(intEstimates));
  relativeErrorTable = zeros(size(intEstimates));

  intTable[:, 1] = intEstimates;

  for k = 2:size(intTable)[1]
    for j = size(intTable)[1]-k+1:-1:1
      intTable[j, k] = ((4^(k-1))intTable[j+1, k-1] - intTable[j, k-1])/((4^(k-1)) - 1);
    end
  end

  for k = 2:size(intTable)[1]
    relativeErrorTable[k, 1] = 100 * abs((intTable[1, k] - intTable[2, k-1]) / intTable[1, k]);
  end

  return (intTable[1, end], relativeErrorTable[1, end]);
end
