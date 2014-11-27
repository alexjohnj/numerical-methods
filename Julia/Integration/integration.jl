function trapezoidInt(f::Function, x1::Number, x2::Number)
  #= Numerical integration of functions using the trapezoid rule

  Keyword arguments:
    f -- Anonymous function to integrate
    x1 -- Lower limit of integration
    x2 -- upper limit of integration

  Returns:
    The the result of ∫f between x1 and x2
  =#
  (x2 - x1) * 0.5 * (f(x1) + f(x2));
end

function trapezoidInt(f::Function, x1::Number, xn::Number, n::Int)
  #= Numerical Integration of functions using the composite trapezoid rule
  Keyword arguments:
    f -- An anonymous function to integrate. Should be able to take vectors
         as arguments.
    x1 -- Lower limit of integration
    xn -- Upper limit of integration
    n -- Number of points to use such that n-1 sections will be used

  Returns:
    The results of ∫f between x1 and x2
  =#
  if n <= 1
    error("Insufficient data points. Try increasing n.");
  end

  h = (xn - x1) / (n - 1);
  0.5 * h * (f(x1) + 2*sum(f([x1+h:h:xn-h])) + f(xn));
end

function simpsonInt(f::Function, x1::Number, x2::Number, x3::Number)
  #= Numerical integration of functions using Simpson's 1/3 Rule
  
  Keyword arguments:
    f -- Function to integrate
    x1 -- Lower limit of integration
    x2 -- Midpoint of integral limits
    x3 -- Upper limit of integration

  Returns:
    The results of ∫f between x1 and x3
  =#
  (x3 - x1) * (1//6) * (f(x1) + 4*f(x2) + f(x3));
end

function simpsonInt(f::Function, x1::Number, xn::Number, n::Int)
  #= Numerical integration using Simpson's Rule
  
  Will select the best method for integration depending on the number of data 
  points provided. If n is odd, Simpson's Composite 1/3 rule will be used. If
  n is a multiple of 3, Simpson's composite 3/8 rule will be used.

  Keyword arguments
    f -- Function to integrate
    x1 -- Lower limit of integration.
    xn -- Upper limit of integration.
    n -- Number of data points to use
  
  Returns
    The results of ∫f between x1 and xn.
  =#
  if n <= 1
    error("Insufficient data points. Try increasing n.");
  elseif isodd(n)
    simpsonOneThirdComposite(f, x1, xn, n);
  elseif n%3 == 0
    simpsonThreeEigthComposite(f, x1, xn, n);
  else
    error("Unusable number of data points. Make sure n is ≥ 1, odd or a multiple of 3.");
  end
end

function simpsonOneThirdComposite(f::Function, x1::Number, xn::Number, n::Int)
  #= Numerical integration using Simpson's 1/3 Rule
  
  Keyword arguments:
    f -- Function to integrate
    x1 -- Lower limit of integration
    xn -- Upper limit of integration
    n -- Number of data points to use. **Must be odd**

  Returns:
    The results of ∫f between x1 and xn
  =#
  if n%2 == 0
    error("Must use an odd number of data points.");
  end 
  h = (xn - x1) / (n - 1);
  (1//3) * h * (f(x1) + 4*(sum(f([x1+h:2h:xn-h]))) + 2*(sum(f([x1+2h:2h:xn-2h]))) + f(xn));
end

function simpsonThreeEigthComposite(f::Function, x1::Number, xn::Number, n::Int)
  #= Numerical integration using Simpson's 3/8 rule

  Keyword arguments:
    f -- Function to integrate
    x1 -- Lower limit of integration
    xn -- Upper limit of integration
    n -- Number of data points to use. **Must be a multiple of 3**.
  
  Returns:
    The result of ∫f between x1 and xn.
  =#
  if n%3 != 0
    error("Must use a multiple of three data points");
  end
  h = (xn - x1) / n;
  (3//8) * h * (f(x1) + 3*sum(f([x1+h:3h:xn-2h])) + 3*sum(f([x1+2h:3h:xn-h])) + 2*sum(f([x1+3h:3h:xn-3h])) + f(xn));
end
