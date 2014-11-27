% Reads data from a plaintext file containing two columns, x and f(x).
% The script then calculates the value f(x) for a user provided number
% by fitting a polynomial of a user provided order.

% Read Data File
fileName = input('Enter filename >> ', 's');
rawData = importdata(fileName);

% Get xint (the number we want to calculate an interpolated value of)
xint = NaN;
while true
    xint = input('Enter the value of x you want to interpolate in GPa >> ');
    if xint > max(rawData(:,1)) || xint < min(rawData(:,1))
       fprintf('%.2f does not lie in the data range\n', xint); 
    else
        break;
    end
end

% Check if f(xint) exists
[containsXint, index] = ismember(xint, rawData(:, 1));
if containsXint
  fprintf('f(%.1f) = %.2fÅ³\n', xint, rawData(index, 2));
  return;
end

% Get the order of the polynomial the user wants to use
polyOrder = input('Enter the order of the polynomial you want to fit >> ');
if polyOrder >= length(rawData)
  error('Insufficient data points. Try reducing the order of the polynomial.')
end

% Find bracketing values
xVals = NaN(polyOrder+1, 1);
fxVals = NaN(polyOrder+1, 1);

for ii = 2:length(rawData)
  if rawData(ii, 1) >= xint
    if ii+floor(polyOrder/2) > length(rawData)
      % If we can't construct symmetric bracketing values around xint
      % because of array bounds errors construct asymmetric ones but 
      % warn the user because I think this is a bad thing but might work.
      leftovers = length(rawData) - (ii + floor(polyOrder/2));
      xVals(:, 1) = rawData(ii-ceil(polyOrder/2)-leftovers:end, 1);
      fxVals(:, 1) = rawData(ii-ceil(polyOrder/2)-leftovers:end, 2);
      warning('Bracketing values are asymmetric');
    elseif ii-ceil(polyOrder/2) < 1
      % Do the same as above but this time in case we try and access
      % negative indices.
      leftovers = 1 - (ii-ceil(polyOrder/2));
      xVals(:, 1) = rawData(1:ii+floor(polyOrder/2)+leftovers, 1);
      fxVals(:, 1) = rawData(1:ii+floor(polyOrder/2)+leftovers, 2);
      warning('Bracketing values are asymmetric');
    else
      % Otherwise construct symmetric bracketing values.
      xVals(:,1) = rawData(ii-ceil(polyOrder/2):ii+floor(polyOrder/2), 1);
      fxVals(:,1) = rawData(ii-ceil(polyOrder/2):ii+floor(polyOrder/2), 2);
    end
    break;
  end
end

% Work out the interpolated value of xint & print it
fxint = newton(xVals, fxVals, xint);
fprintf('f(%.1f) = %.2fÅ³\n', xint, fxint);
