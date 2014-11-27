function fixedpointiteration(g::Function, xi::Number, ea::Number, maxIter::Int64)
  # Estiamte roots for an equation using the fixed-point iteration method
  for i = 1:maxIter
    xiold = xi;
    xi = g(xi);

    if (abs(xi - xiold) / xi) < ea
      break;
    end;
  end
  
  return xi;
end
