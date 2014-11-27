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
