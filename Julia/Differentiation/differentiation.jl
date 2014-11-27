function ffdiff(f::Function, x1::Number, h::Number)
  # Estimates the derivative of a function using a forward finite difference
  # approximation
  (f(x1 + h) - f(x1)) / (x1+h - x1);
end

function bfdiff(f::Function, x1::Number, h::Number)
  # Estimates the derivative of a function using a backwards finite difference
  # approximation
  (f(x1) - f(x1-h)) / (x1 + h - x1);
end

function cfdiff(f::Function, x1::Number, h::Number)
  # Estimates the derivative of a function using a centred finite difference
  # approximation  
  (f(x1 + h)  - f(x1 - h)) / 2(x1 + h - x1);
end