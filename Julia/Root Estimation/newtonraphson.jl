function newtonRaphson(f::Function, fprime::Function, xi::Number, ea::Number, maxIter::Int64)
  # Estimate roots of functions using the Newton-Raphson method
  for i = 1:maxIter
    xiprev  = xi;
    xi = xi - (f(xi) / fprime(xi));

    if abs(xi - xiprev) / xi < ea
      break;
    end;
  end
  return xi;
end
