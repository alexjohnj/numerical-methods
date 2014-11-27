function golden(f::Function, xl::Number, xu::Number, ea::Number, maxIter::Int64)
  # Estimate the optimum of a function using the golden-section search method
  # Interesting note, removing array allocations from this function improved execution
  # time by a factor of 100!

  const ϕ = (1+sqrt(5)) / 2; # Hope you like unicode 'cause I'm using it
  xo = NaN; # Current estimate for the optimum

  # Calculate x1, x2 & f(x1), f(x2) for the initial lower and upper bounds
  x1 = xl + (xu - xl) / ϕ; # Upper test point
  x2 = xu - (xu - xl) / ϕ; # Lower test point

  f1 = f(x1);
  f2 = f(x2);

  for i = 1:maxIter
    if f1 > f2
      # The optimum is to the left of x1. Make the upper bound x1, set the
      # upper test point to the previous lower test point and calculate a
      # new lower test point for the new upper bound.
      xo = x2;
      xu = x1;
      x1 = x2;
      f1 = f2;
      x2 = xu - (xu - xl) / ϕ;
      f2 = f(x2);

      err = abs((x2 - xl) / xo) * 100;

      if err < ea
        break;
      end

    else
      # The optimum is to the right of x1. Make the lower bound x2, set the
      # lower test point to the previous upper test point and calculate a
      # new upper test point for the new lower bounds.
      xo = x1;
      xl = x2;
      x2 = x1;
      f2 = f1;
      x1 = xl + (xu - xl) / ϕ;
      f1 = f(x1);

      err = abs((x2 - xl) / xo) * 100;

      if err < ea
        break;
      end
    end
  end
  return xo;
end
